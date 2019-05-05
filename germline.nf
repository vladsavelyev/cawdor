#!/usr/bin/env nextflow

if (params.help) exit 0, helpMessage()

genomeFasta = Utils.findRefFile(params, 'genomeFasta')
genomeIndex = Utils.findRefFile(params, 'genomeIndex')
genomeDict  = Utils.findRefFile(params, 'genomeDict')
dbsnp       = Utils.findRefFile(params, 'dbsnp')
dbsnpIndex  = Utils.findRefFile(params, 'dbsnpIndex')
intervals   = Utils.findRefFile(params, 'intervals')

// Check for awsbatch profile configuration
// make sure queue is defined
if (workflow.profile == 'awsbatch') {
    if (!params.awsqueue) exit 1, "Provide the job queue for aws batch!"
}

// Set up the bamFiles channel
if (params.containsKey("samples")) {
  tsvFile = file(params.samples)
  if (!file(tsvFile).exists()) {
    exit 1, "Input --samples file ${tsvFile} does not exist"
  }
} else {
  tsvFile = "${params.outDir}/Align/samples.tsv"
  if (!file(tsvFile).exists()) {
    exit 1, "${tsvFile} not found in output directory. Make sure you ran align.nf on this folder. See --help"
  }
}
inputData = Utils.extractBams(file(tsvFile))
inputChannel = Channel.from(inputData)
// We need only normal BAM
inputChannel = inputChannel.map {
  idPatient, idSampleNormal, bamNormal, baiNormal, idSampleTumour, bamTumour, baiTumour ->
    [idSampleNormal, bamNormal, baiNormal]
}

Utils.startMessage(log, workflow, config, params, tsvFile)

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

// We know that MuTect2 (and other somatic callers) are notoriously slow.
// To speed them up we are chopping the reference into smaller pieces.
// Do variant calling by this intervals, and re-merge the VCFs.
// Since we are on a cluster or a multi-CPU machine, this can parallelize the
// variant call processes and push down the variant call wall clock time significanlty.
process CreateIntervalBeds {
  tag {intervals.fileName}

  input:
  file(intervals) from Channel.value(intervals)

  output:
  file '*.bed' into intervalBeds mode flatten

  script:
  // If the interval file is BED format, the fifth column is interpreted to
  // contain runtime estimates, which is then used to combine short-running jobs
  if (intervals.getName().endsWith('.bed'))
    """
    awk -v FS="\t" '{
      t = \$5  # runtime estimate
      if (t == "") {
        # no runtime estimate in this row, assume default value
        t = (\$3 - \$2) / ${params.nucleotidesPerSecond}
      }
      if (chunk_fname == "" || (chunk > 600 && (chunk + t) > longest * 1.05)) {
        # start a new chunk
        chunk_fname = sprintf("%s_%d-%d.bed", \$1, \$2+1, \$3)
        chunk = 0
        longest = 0
      }
      if (t > longest)
        longest = t
      chunk += t
      print \$0 > chunk_fname
      close(chunk_fname)
    }' ${intervals}
    """
  else
    """
    awk -v FS="[:-]" '{
      chunk_fname = sprintf("%s_%d-%d.bed", \$1, \$2, \$3);
      printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > chunk_fname
      close(chunk_fname)
    }' ${intervals}
    """
}

intervalBeds = intervalBeds.ifEmpty{ exit 1, "No intervals" }

intervalBeds = intervalBeds
  .map { intervalBed ->
    def duration = 0.0
    for (line in intervalBed.readLines()) {
      final fields = line.split('\t')
      if (fields.size() >= 5) duration += fields[4].toFloat()
      else {
        start = fields[1].toInteger()
        end = fields[2].toInteger()
        duration += (end - start) / params.nucleotidesPerSecond
      }
    }
    [duration, intervalBed]
  }.toSortedList({ a, b -> b[0] <=> a[0] })
    .flatten().collate(2)
    .map{duration, intervalFile -> intervalFile}

if (params.verbose) bedIntervals = bedIntervals.view {
  "  Interv: ${it.baseName}"
}

if (params.verbose) inputChannel = inputChannel.view {
  "Normal sample: ${it[0]}, BAM: [${it[1]}"
}

// Manta, Strelka (no splitting by intervals)
inputChannel.into { bamsForManta; bamsForStrelka; bamsForIntervals }

// MuTect, VarDict (partitioning by intervals)
bamsByIntervals = bamsForIntervals.combine(intervalBeds)  // cartesian product of 2 channels
bamsByIntervals.into { bamsForHC; bamsForVardict }

process RunHaplotypecaller {
  tag {idSample + "-haplotypecaller-" + intervalBed.baseName}

  input:
  set idSample, file(bam), file(bai), file(intervalBed) from bamsForHC
  set file(genomeFasta), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex) from Channel.value([
    genomeFasta,
    genomeIndex,
    genomeDict,
    dbsnp,
    dbsnpIndex
  ])

  output:
  set val("haplotypecaller"), idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into hcOutput

  when: !params.strelkaOnly && !params.onlyQC

  script:
  outFile = "${idSample}-interval_${intervalBed.baseName}.vcf.gz"
  """
  gatk --java-options "-Xms6g -Xmx${task.memory.toGiga()}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
  HaplotypeCaller \
  -R ${genomeFasta} \
  -I ${bam} \
  -L ${intervalBed} \
  -D ${dbsnp} \
  -O ${outputFile} \
  --annotation MappingQualityRankSumTest \
  --annotation MappingQualityZero \
  --annotation QualByDepth \
  --annotation ReadPosRankSumTest \
  --annotation RMSMappingQuality \
  --annotation BaseQualityRankSumTest \
  --annotation FisherStrand \
  --annotation MappingQuality \
  --annotation DepthPerAlleleBySample \
  --annotation Coverage \
  --annotation ClippingRankSumTest \
  --annotation DepthPerSampleHC \
  --interval-set-rule INTERSECTION \
  --native-pair-hmm-threads 1 \
  """
}
hcGenomicVCF = hcGenomicVCF.groupTuple(by:[0,1,2,3])

process RunVarDictGermline {
  tag {idSample + "-vardict-germline" + intervalBed.baseName}

  input:
  set idSample, file(bam), file(bai), file(intervalBed) from bamsForVardict
  set file(genomeFasta), file(genomeIndex), file(genomeDict) from Channel.value([
      genomeFasta,
      genomeIndex,
      genomeDict
  ])

  output:
  set val("vardict"), idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into vardictOutput

  when: !params.strelkaOnly && !params.onlyQC

  script:
  tmpDir = file("tmp").mkdir()
  outFile = "${idSample}-interval_${intervalBed.baseName}.vcf.gz"
  """
  unset JAVA_HOME && \
  export VAR_DICT_OPTS='-Xms750m -Xmx${task.memory.toGiga()}g -XX:+UseSerialGC -Djava.io.tmpdir=${tmpDir}' && \
  vardict-java -G ${genomeFasta} \
  -N ${idSample} \
  -b ${bam} \
  -c 1 -S 2 -E 3 -g 4 --nosv --deldupvar -Q 10 -F 0x700 -f 0.1 \
  ${intervalBed} \
  | teststrandbias.R \
  | var2vcf_valid.pl -A -N ${idSample} -E -f 0.1
  | bcftools filter -i 'QUAL >= 0' \
  | bcftools filter --soft-filter 'LowFreqBias' --mode '+' -e 'FORMAT/AF[0] < 0.02 && FORMAT/VD[0] < 30 \
  && INFO/SBF < 0.1 && INFO/NM >= 2.0' \
  | awk -F\$'\t' -v OFS='\t' '{if (\$0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", \$4) } {print}' \
  | awk -F\$'\t' -v OFS='\t' '{if (\$0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", \$5) } {print}' \
  | awk -F\$'\t' -v OFS='\t' '\$1!~/^#/ && \$4 == \$5 {next} {print}' \
  | vcfstreamsort \
  | bgzip -c > ${outFile} && tabix -p vcf ${outFile} \
  """
}

vardictOutput = vardictOutput.groupTuple(by:[0,1,2,3])

// we are merging the VCFs that are called separately for different intervals
// so we can have a single sorted VCF containing all the calls for a given caller
vcfsToConcat = hcOutput.mix(vardictOutput)
if (params.verbose) vcfsToConcat = vcfsToConcat.view {
  "Interval VCFs to be concatenated:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  Files : ${it[4].fileName}"
}
process ConcatVCFbyInterval {
  tag {idSample + "-" + variantCaller}

  publishDir "${params.outDir}/VariantCalling/${idSample}/${"$variantCaller"}", mode: params.publishDirMode

  input:
  set variantCaller, idSample, file(vcfFiles), file(tbiFiles) from vcfsToConcat

  output:
  set variantCaller, idSample,
      file("*-${variantCaller}.vcf.gz"), file("*-${variantCaller}.vcf.gz.tbi") into vcfConcatenated

  when: !params.onlyQC

  script:
  outFile = "${idSample}-germline-${variantCaller}.vcf.gz"  // will add .gz and .gz.tbi automatically
  """
  bcftools concat ${idSample}-interval_*.vcf.gz -a | bcftools sort -Oz -o ${outFile}
  tabix -p vcf ${outFile}
  """
}

if (params.verbose) vcfConcatenated = vcfConcatenated.view {
  "Variant calling output:\n\
  Tool: ${it[0]}\tSample: ${it[1]}\tVCF : ${it[2].fileName}"
}

process RunManta {
  tag {idSample}

  publishDir "${params.outDir}/VariantCalling/${idSample}/Manta", mode: params.publishDirMode

  input:
  set idSample, file(bam), file(bai) from bamsForManta
  file(targetBED) from Channel.value(params.targetBED ? file(params.targetBED) : "null")
  set file(genomeFasta), file(genomeIndex) from Channel.value([
      genomeFasta,
      genomeIndex
  ])

  output:
  set val("manta"), idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into mantaOutput

  when: !params.onlyQC

  script:
  targetCmd = params.targetBED ?
      "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" :
      ""
  options = params.targetBED ? "--exome --callRegions call_targets.bed.gz" : ""
  outputFile = "${idSample}-germline-manta.vcf.gz"
  """
  ${targetCmd}
  configManta.py \
  --bam ${bam} \
  --reference ${genomeFasta} \
  ${options} \
  --runDir Manta

  python Manta/runWorkflow.py -m local -j ${task.cpus}

  mv Manta/results/variants/diploidSV.vcf.gz ${outputFile}
  mv Manta/results/variants/diploidSV.vcf.gz.tbi ${outputFile}.tbi
  """
}

if (params.verbose) mantaOutput = mantaOutput.view {
  "Variant calling output:\n\
  Tool: ${it[0]}\tSample: ${it[1]}\tVCF : ${it[2].fileName}"
}

process RunStrelka {
  tag {idSample}

  publishDir "${params.outDir}/VariantCalling/${idSample}/Strelka", mode: params.publishDirMode

  input:
  set idSample, file(bam), file(bai) from bamsForStrelka
  file(targetBED) from Channel.value(params.targetBED ? file(params.targetBED) : "null")
  set file(genomeFasta), file(genomeIndex) from Channel.value([
      genomeFasta,
      genomeIndex
  ])

  output:
  set val("strelka"), idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into strelkaOutput

  when: !params.onlyQC

  script:
  targetsCmd = params.targetBED ?
      "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" :
      ""
  options = params.targetBED ? "--exome --callRegions call_targets.bed.gz" : ""
  outputFile = "${idSample}-germline-strelka.vcf.gz"
  """
  ${targetsCmd}
  configureStrelkaGermlineWorkflow.py \
  --bam ${bam} \
  --referenceFasta ${genomeFasta} \
  ${options} \
  --runDir Strelka

  python Strelka/runWorkflow.py -m local -j ${task.cpus}
  mv Strelka/results/variants/variants.vcf.gz ${outputFile}
  mv Strelka/results/variants/variants.vcf.gz.tbi ${outputFile}.tbi
  """
}

if (params.verbose) strelkaOutput = strelkaOutput.view {
  "Variant calling output:\n\
  Tool: ${it[0]}\tSample: ${it[1]}\tVCF : ${it[2].fileName}"
}

vcfForQC = Channel.empty().mix(
  vcfConcatenated.map {
    variantcaller, idSample, vcf, tbi ->
    [variantcaller, vcf]
  },
  strelkaOutput.map {
    variantcaller, idSample, vcf, tbi ->
    [variantcaller, vcf[1]]
  },
  mantaOutput.map {
    variantcaller, idSample, vcf, tbi ->
    [variantcaller, vcf[2]]
  })

(vcfForBCFtools, vcfForVCFtools) = vcfForQC.into(2)

process RunBcftoolsStats {
  tag {vcf}

  publishDir "${params.outDir}/Reports/BCFToolsStats", mode: params.publishDirMode

  input:
    set variantCaller, file(vcf) from vcfForBCFtools

  output:
  file ("*_stats.txt") into bcfReport

  script:
  """
  bcftools stats ${vcf} > ${vcf.simpleName}.bcftools_stats.txt
  """
}

if (params.verbose) bcfReport = bcfReport.view {
  "BCFTools stats report:\n\
  File  : [${it.fileName}]"
}

process RunVcftools {
  tag {vcf}

  publishDir "${params.outDir}/Reports/VCFTools", mode: params.publishDirMode

  input:
    set variantCaller, file(vcf) from vcfForVCFtools

  output:
    file ("${vcf.simpleName}.*") into vcfReport

  script:
  """
  vcftools \
  --gzvcf ${vcf} \
  --relatedness2 \
  --out ${vcf.simpleName}

  vcftools \
  --gzvcf ${vcf} \
  --TsTv-by-count \
  --out ${vcf.simpleName}

  vcftools \
  --gzvcf ${vcf} \
  --TsTv-by-qual \
  --out ${vcf.simpleName}

  vcftools \
  --gzvcf ${vcf} \
  --FILTER-summary \
  --out ${vcf.simpleName}
  """
}

if (params.verbose) vcfReport = vcfReport.view {
  "VCFTools stats report:\n\
  File  : [${it.fileName}]"
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def helpMessage() {
  // Display help message
  Utils.cawdorMsg(log)
  log.info "    Usage:"
  log.info "       nextflow run germline.nf --samples <file.tsv> --outDir <Dir> --genome <Genome>"
  log.info ""
  log.info "    --outDir <Directory>"
  log.info "       Output directory. Can be the output from align.nf subworkflow, "
  log.info "       in which case it will use Align/sample.tsv if --samples not specified"
  log.info "    --samples <file.tsv>"
  log.info "       Specify a TSV file containing paths to sample files."
  log.info "    --test"
  log.info "       Use a test sample."
  log.info "    --genome <Genome>"
  log.info "       Use a specific genome version."
  log.info "       Possible values are:"
  log.info "         GRCh37"
  log.info "         GRCh38 (Default)"
  log.info "         smallGRCh37 (Use a small reference (Tests only))"
  log.info "    --onlyQC"
  log.info "       Run only QC tools and gather reports"
  log.info "    --help"
  log.info "       you're reading it"
  log.info "    --verbose"
  log.info "       Adds more verbosity to workflow"
}

workflow.onComplete {
  Utils.endMessage(log, workflow, config, params, tsvFile)
}

workflow.onError {
  // Display error message
  Utils.nextflowMessage(log, workflow)
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
