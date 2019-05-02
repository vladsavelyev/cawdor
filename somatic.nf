#!/usr/bin/env nextflow

if (params.help) exit 0, helpMessage()

genomeFasta = Utils.findRefFile(params, 'genomeFasta')
genomeIndex = Utils.findRefFile(params, 'genomeIndex')
genomeDict  = Utils.findRefFile(params, 'genomeDict')
purpleHet   = Utils.findRefFile(params, 'purpleHet')
purpleGC    = Utils.findRefFile(params, 'purpleGC')
intervals   = Utils.findRefFile(params, 'intervals')
callable    = Utils.findRefFile(params, 'callable')
if (![genomeFasta, genomeIndex, genomeDict, intervals].every())
  exit 1, "Missing reference files for alignment. See --help for more information"

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

Utils.startMessage(log, workflow, config, params)
log.info "Callable    : " + callable

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

if (params.verbose) intervalBeds = intervalBeds.view {
  "  Interv: ${it.baseName}"
}

if (params.verbose) inputChannel = inputChannel.view {
  "BAMs for variant calling:\n\
  ID    : ${it[0]},\tSamples: [${it[1]}, ${it[4]}]\n\
  Files : [${it[2].fileName}, ${it[5].fileName}]"
}

// Manta, Strelka (no splitting by intervals)
inputChannel.into { bamsForManta; bamsForStrelkaBP; bamsForCobalt; bamsForAmber; bamsForIntervals }

// MuTect, VarDict (partitioning by intervals)
bamsByIntervals = bamsForIntervals.combine(intervalBeds)  // cartesian product of 2 channels
bamsByIntervals.into { bamsForMutect; bamsForVardict }

// This will give as a list of unfiltered calls for Mutect.
process RunMutect {
  tag {idPatient + "-mutect-" + intervalBed.baseName}

  input:
  set idPatient, idSampleNormal, file(bamNormal), file(baiNormal),
      idSampleTumour, file(bamTumour), file(baiTumour), file(intervalBed) from bamsForMutect
  set file(genomeFasta), file(genomeIndex), file(genomeDict) from Channel.value([
      genomeFasta, genomeIndex, genomeDict])
  file(callable) from Channel.value(callable ? file(callable) : "null")

  output:
  set val("mutect"), idPatient, idSampleNormal, idSampleTumour,
      file("*.vcf.gz"), file("*.vcf.gz.tbi") into mutectOutput

  when: !params.strelkaOnly && !params.onlyQC

  script:
  def makeTarget = ""
  def targetBed = intervalBed
  if (callable) {
    targetBed = 'target.bed'
    makeTarget = "bedtools intersect -a ${intervalBed} -b ${callable} > ${targetBed}"
  }
  outFile = "${idPatient}-interval_${intervalBed.baseName}.vcf.gz"
  """
  ${makeTarget}
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
  Mutect2 \
  -R ${genomeFasta}\
  -I ${bamTumour}  -tumor ${idSampleTumour} \
  -I ${bamNormal} -normal ${idSampleNormal} \
  -L ${targetBed} \
  -O ${outFile} \
  """
}

mutectOutput = mutectOutput.groupTuple(by:[0,1,2,3])  // group by normals

//process RunFreeBayes {
//  tag {idSampleTumour + "_vs_" + idSampleNormal + "-" + intervalBed.baseName}
//
//  input:
//    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal),
//                   idSampleTumour, file(bamTumour), file(baiTumour), file(intervalBed) from bamsForFreebayes
//    file(genomeFasta) from Channel.value(genomeFasta)
//    file(genomeIndex) from Channel.value(genomeIndex)
//
//  output:
//    set val("FreeBayes"), idPatient, idSampleNormal, idSampleTumour,
//        file("${intervalBed.baseName}_${idSampleTumour}_vs_${idSampleNormal}.vcf") into freebayesOutput
//
//  when: !params.onlyQC
//
//  script:
//  """ \
//  freebayes \
//  -f ${genomeFasta} \
//  --pooled-continuous \
//  --pooled-discrete \
//  --genotype-qualities \
//  --report-genotype-likelihood-max \
//  --allele-balance-priors-off \
//  --min-alternate-fraction 0.03 \
//  --min-repeat-entropy 1 \
//  --min-alternate-count 2 \
//  -t ${intervalBed} \
//  ${bamTumour} \
//  ${bamNormal} > ${intervalBed.baseName}_${idSampleTumour}_vs_${idSampleNormal}.vcf \
//  """
//}
//
//freebayesOutput = freebayesOutput.groupTuple(by:[0,1,2,3])

process RunVarDict {
  tag {idPatient + "-vardict-" + intervalBed.baseName}

  input:
  set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumour, file(bamTumour),
    file(baiTumour), file(intervalBed) from bamsForVardict
  file(genomeFasta) from Channel.value(genomeFasta)
  file(genomeIndex) from Channel.value(genomeIndex)
  file(callable) from Channel.value(callable ? file(callable) : "null")

  output:
  set val("vardict"), idPatient, idSampleNormal, idSampleTumour,
      file("*.vcf.gz"), file("*.vcf.gz.tbi") into vardictOutput

  when: !params.strelkaOnly && !params.onlyQC

  script:
  def makeTarget = ""
  def targetBed = intervalBed
  if (callable) {
    targetBed = 'target.bed'
    makeTarget = "bedtools intersect -a ${intervalBed} -b ${callable} > ${targetBed}"
  }
  tmpDir = file("tmp").mkdir()
  outFile = "${idPatient}-interval_${intervalBed.baseName}.vcf.gz"
  """
  ${makeTarget}
  unset JAVA_HOME && \
  export VAR_DICT_OPTS='-Xms750m -Xmx${task.memory.toGiga()}g -XX:+UseSerialGC -Djava.io.tmpdir=${tmpDir}' && \
  vardict-java -G ${genomeFasta} \
  -N ${idSampleNormal} \
  -b "${bamTumour}|${bamNormal}" \
  -c 1 -S 2 -E 3 -g 4 --nosv --deldupvar -Q 10 -F 0x700 -f 0.1 \
  ${targetBed} \
  | awk 'NF>=48' \
  | testsomatic.R \
  | var2vcf_paired.pl -P 0.9 -m 4.25 -f 0.01 -M  -N "${idSampleTumour}|${idSampleNormal}" \
  | bcftools filter -m '+' -s 'REJECT' -e 'STATUS !~ ".*Somatic"' \
  2> /dev/null \
  | bcftools filter --soft-filter 'LowFreqBias' --mode '+' -e 'FORMAT/AF[0] < 0.02 && \
  FORMAT/VD[0] < 30 && FORMAT/SBF[0] < 0.1 && FORMAT/NM[0] >= 2.0' \
  | bcftools filter -i 'QUAL >= 0' \
  | sed 's/\\\\.*Somatic\\\\/Somatic/' \
  | sed 's/REJECT,Description=".*">/REJECT,Description="Not Somatic via VarDict">/' \
  | awk -F\$'\\t' -v OFS='\\t' '{if (\$0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", \$4) } {print}' \
  | awk -F\$'\\t' -v OFS='\\t' '{if (\$0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", \$5) } {print}' \
  | awk -F\$'\\t' -v OFS='\\t' '\$1!~/^#/ && \$4 == \$5 {next} {print}' \
  | vcfstreamsort \
  | bgzip -c > ${outFile} && tabix -p vcf ${outFile} \
  """
}

vardictOutput = vardictOutput.groupTuple(by:[0,1,2,3])

// we are merging the VCFs that are called separately for different intervals
// so we can have a single sorted VCF containing all the calls for a given caller
vcfsToConcat = mutectOutput.mix(vardictOutput)
if (params.verbose) vcfsToConcat = vcfsToConcat.view {
  "Interval VCFs to be concatenated:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  Files : ${it[4].fileName}"
}
process ConcatVCFbyInterval {
  tag {idPatient + "-" + variantCaller}

  publishDir "${params.outDir}/VariantCalling/${idPatient}/${"$variantCaller"}", mode: params.publishDirMode

  input:
  set variantCaller, idPatient, idSampleNormal, idSampleTumour,
      file(vcfFiles), file(tbiFiles) from vcfsToConcat

  output:
  set variantCaller, idPatient, idSampleNormal, idSampleTumour,
      file("*-${variantCaller}.vcf.gz"), file("*-${variantCaller}.vcf.gz.tbi") into vcfConcatenated

  when: !params.onlyQC

  script:
  outFile = "${idPatient}-${variantCaller}.vcf.gz"  // will add .gz and .gz.tbi automatically
  """
  bcftools concat ${idPatient}-interval_*.vcf.gz -a | bcftools sort -Oz -o ${outFile}
  tabix -p vcf ${outFile}
  """
}

if (params.verbose) vcfConcatenated = vcfConcatenated.view {
  "Variant calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  File  : ${it[4].fileName}"
}

process RunManta {
  tag {idSampleTumour + "_vs_" + idSampleNormal}

  publishDir "${params.outDir}/VariantCalling/${idPatient}/Manta", mode: params.publishDirMode

  input:
  set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumour, file(bamTumour),
    file(baiTumour) from bamsForManta
  set file(genomeFasta), file(genomeIndex) from Channel.value([genomeFasta, genomeIndex])

  output:
  set val("manta"), idPatient, idSampleNormal, idSampleTumour, file("*.vcf.gz"), file("*.vcf.gz.tbi") into mantaOutput
  set idPatient, idSampleNormal, idSampleTumour, \
      file("Manta/results/variants/candidateSmallIndels.vcf.gz"), \
      file("Manta/results/variants/candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

  when: !params.onlyQC

  script:
  outputFile = "${idPatient}-manta.vcf.gz"
  """ \
  configManta.py \
  --normalBam ${bamNormal} \
  --tumorBam ${bamTumour} \
  --reference ${genomeFasta} \
  --runDir Manta

  python Manta/runWorkflow.py -m local -j ${task.cpus}

  mv Manta/results/variants/somaticSV.vcf.gz ${outputFile}
  mv Manta/results/variants/somaticSV.vcf.gz.tbi ${outputFile}.tbi
  """
}

if (params.verbose) mantaOutput = mantaOutput.view {
  "Manta output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  Files : ${it[4].fileName}"
}

// Running Strelka Best Practice with Manta indel candidates
// For easier joining, remaping channels to idPatient, idSampleNormal, idSampleTumour...
bamsForStrelkaBP = bamsForStrelkaBP.map {
   idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumour, bamTumour, baiTumour ->
  [idPatientNormal, idSampleNormal, idSampleTumour, bamNormal, baiNormal, bamTumour, baiTumour]
}.join(mantaToStrelka, by:[0,1,2]).map {
   idPatientNormal, idSampleNormal, idSampleTumour, bamNormal, baiNormal, bamTumour, baiTumour, mantaCSI, mantaCSIi ->
  [idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumour, bamTumour, baiTumour, mantaCSI, mantaCSIi]
}

process RunStrelka {
  tag {idSampleTumour + "_vs_" + idSampleNormal}

  publishDir "${params.outDir}/VariantCalling/${idPatient}/Strelka", mode: params.publishDirMode

  input:
  set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumour, file(bamTumour), file(baiTumour),
      file(mantaCSI), file(mantaCSIi) from bamsForStrelkaBP
  set file(genomeFasta), file(genomeIndex), file(genomeDict) from Channel.value([
    genomeFasta,
    genomeIndex,
    genomeDict
  ])

  output:
  set val("Strelka"), idPatient, idSampleNormal, idSampleTumour,
      file("*.vcf.gz"), file("*.vcf.gz.tbi") into strelkaOutput

  when: !params.onlyQC

  script:
  outputFile = "${idPatient}-strelka.vcf.gz"
  """ 
  configureStrelkaSomaticWorkflow.py \
  --tumor ${bamTumour} \
  --normal ${bamNormal} \
  --referenceFasta ${genomeFasta} \
  --indelCandidates ${mantaCSI} \
  --runDir Strelka

  python Strelka/runWorkflow.py -m local -j ${task.cpus}

  bcftools concat -a Strelka/results/variants/somatic.indels.vcf.gz Strelka/results/variants/somatic.snvs.vcf.gz \
  -Oz -o ${outputFile}
  tabix -p vcf ${outputFile}
  """
}

if (params.verbose) strelkaOutput = strelkaOutput.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  File  : ${it[4].fileName}"
}

process RunAmber {
  label "purple"

  tag {idSampleTumour + "_vs_" + idSampleNormal}

  publishDir "${params.outDir}/VariantCalling/${idPatient}/Purple/amber", mode: params.publishDirMode

  input:
  set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumour, file(bamTumour), file(baiTumour) from bamsForAmber
  set file(genomeFasta), file(genomeIndex), file(genomeDict), file(purpleHet) from Channel.value([
    genomeFasta,
    genomeIndex,
    genomeDict,
    purpleHet
  ])

  output:
  set idPatient, idSampleNormal, idSampleTumour, file('amber/*') into amberOutput

  when: !params.onlyQC

  script:
  jvm_opts = "-Xms750m -Xmx${task.memory.toGiga()}g"
  """ \
  AMBER ${jvm_opts} \
  -tumor ${idPatient} \
  -tumor_bam ${bamTumour} \
  -reference ${idSampleNormal} \
  -reference_bam ${bamNormal} \
  -ref_genome ${genomeFasta} \
  -bed ${purpleHet} \
  -threads ${task.cpus} \
  -output_dir amber \
  """
}

process RunCobalt {
  label "purple"

  tag {idSampleTumour + "_vs_" + idSampleNormal}

  publishDir "${params.outDir}/VariantCalling/${idPatient}/Purple/cobalt", mode: params.publishDirMode

  input:
  set idPatient, idSampleNormal, file(bamNormal), file(baiNormal),
      idSampleTumour, file(bamTumour), file(baiTumour) from bamsForCobalt
  set file(genomeFasta), file(genomeIndex), file(genomeDict), file(purpleGC) from Channel.value([
      genomeFasta,
      genomeIndex,
      genomeDict,
      purpleGC
  ])

  output:
  set idPatient, idSampleNormal, idSampleTumour, file('cobalt/*') into cobaltOutput

  when: !params.onlyQC

  script:
  jvm_opts = "-Xms750m -Xmx${task.memory.toGiga()}g"
  """ \
  COBALT ${jvm_opts} \
  -tumor ${idPatient} \
  -tumor_bam ${bamTumour} \
  -reference ${idSampleNormal} \
  -reference_bam ${bamNormal} \
  -ref_genome ${genomeFasta} \
  -gc_profile ${purpleGC} \
  -threads ${task.cpus} \
  -output_dir cobalt \
  """
}

mantaOutput.into { mantaOutputForQC; mantaOutputForPurple }

process RunPurple {
  label "purple"

  tag {idSampleTumour + "_vs_" + idSampleNormal}

  publishDir "${params.outDir}/VariantCalling/${idPatient}/Purple/purple", mode: params.publishDirMode

  input:
  set idPatient, idSampleNormal, idSampleTumour, file(amberFiles) from amberOutput
  set idPatient, idSampleNormal, idSampleTumour, file(cobaltFiles) from cobaltOutput
  set caller, idPatient, idSampleNormal, idSampleTumour, file(mantaVCF), file(mantaIndex) from mantaOutputForPurple
  set file(genomeFasta), file(genomeIndex), file(genomeDict), file(purpleGC) from Channel.value([
      genomeFasta,
      genomeIndex,
      genomeDict,
      purpleGC
  ])

  output:
  set val("PURPLE"), idPatient, idSampleNormal, idSampleTumour, file('purple/*') into purpleOutput
//      '*.purple.cnv', '*.purple.gene.cnv',
//      '*.purple.germline.cnv', '*.purple.sv.vcf.gz', '*.purple.sv.vcf.gz.tbi', '*.purple.purity',
//      '*.purple.qc' into purpleOutput

  when: !params.onlyQC

  script:
  jvm_opts = "-Xms750m -Xmx${task.memory.toGiga()}g"
//  -somatic_vcf ${somaticVCF}
//  - circos circos
  """ \
  mkdir -p purple/amber purple/cobalt
  mv ${amberFiles} purple/amber
  mv ${cobaltFiles} purple/cobalt
  PURPLE ${jvm_opts} \
  -run_dir purple \
  -output_dir purple \
  -reference ${idSampleNormal} \
  -tumor ${idPatient} \
  -threads ${task.cpus} \
  -gc_profile ${purpleGC} \
  -structural_vcf ${mantaVCF} \
  """
}

def vcfForQC = Channel.empty().mix(
  vcfConcatenated.map {
    variantcaller, idPatient, idSampleNormal, idSampleTumour, vcf, tbi ->
    [variantcaller, vcf]
  },
  mantaOutputForQC.map {
    variantcaller, idPatient, idSampleNormal, idSampleTumour, vcf, tbi ->
    [variantcaller, vcf]
  },
  strelkaOutput.map {
    variantcaller, idPatient, idSampleNormal, idSampleTumour, vcf, tbi ->
    [variantcaller, vcf]
  }
)

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
  log.info "       nextflow run somatic.nf --samples <file.tsv> --outDir <Dir> --genome <Genome>"
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
  log.info "    --strelkaOnly"
  log.info "       Run only Strelka variant caller (to avoid extensive run times for FFPE samples)"
  log.info "    --help"
  log.info "       you're reading it"
  log.info "    --verbose"
  log.info "       Adds more verbosity to workflow"
}

workflow.onComplete {
  Utils.endMessage(log, workflow, config, params)
}

workflow.onError {
  // Display error message
  Utils.nextflowMessage(log, workflow)
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
