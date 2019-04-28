#!/usr/bin/env nextflow

if (params.help) exit 0, helpMessage()

genomeFasta = Utils.findRefFile(params, 'genomeFasta')
genomeIndex = Utils.findRefFile(params, 'genomeIndex')
genomeDict  = Utils.findRefFile(params, 'genomeDict')
dbsnp       = Utils.findRefFile(params, 'dbsnp')
intervals   = Utils.findRefFile(params, 'intervals')
if (![genomeFasta, genomeIndex, genomeDict, dbsnp, intervals].every())
  exit 1, "Missing reference files for alignment. See --help for more information"

// Set up the bamFiles channel
bamFiles = Channel.empty()
if (params.containsKey("samples")) {
  tsvFile = file(params.samples)
  if (!file(tsvFile).exists()) {
    exit 1, "Input --samples file ${tsvFile} does not exist"
  }
} else {
  tsvFile = "${params.outDir}/samples.tsv"
  if (!file(tsvFile).exists()) {
    exit 1, "${tsvFile} not found in output directory. Make sure you ran align.nf on this folder. See --help"
  }
}
bamFiles = Utils.extractBams(file(tsvFile))

minimalInformationMessage()

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

if (params.verbose) bamFiles = bamFiles.view {
  "BAMs for variant calling:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  Files : [${it[3].fileName}, ${it[4].fileName}]"
}

// separate recalibrateBams by status
bamsNormal = Channel.create()
bamsTumour = Channel.create()

bamFiles.choice(bamsTumour, bamsNormal) {it[1] == 0 ? 1 : 0}

bamsNormal = bamsNormal.ifEmpty{exit 1, "No normal sample defined, check TSV file: ${tsvFile}"}
bamsTumour = bamsTumour.ifEmpty{exit 1, "No tumour sample defined, check TSV file: ${tsvFile}"}

// Removing status because not relevant anymore
bamsNormal = bamsNormal.map { idPatient, status, idSample, bam, bai -> [idPatient, idSample, bam, bai] }
bamsTumour = bamsTumour.map { idPatient, status, idSample, bam, bai -> [idPatient, idSample, bam, bai] }

// We know that MuTect2 (and other somatic callers) are notoriously slow.
// To speed them up we are chopping the reference into smaller pieces.
// Do variant calling by this intervals, and re-merge the VCFs.
// Since we are on a cluster or a multi-CPU machine, this can parallelize the
// variant call processes and push down the variant call wall clock time significanlty.

params.str = 'Hello world!'

//process splitLetters {
//
//  output:
//  file 'chunk_*' into letters mode flatten
//
//  """
//  printf '${params.str}' | split -b 6 - chunk_
//  """
//}
//
//if (params.verbose) letters = letters.view {
//  "  Letter: ${it.baseName}"
//}

process CreateIntervalBeds {
  tag {intervals.fileName}

  input:
    file(intervals) from Channel.value(intervals)

  output:
    file '*.bed' into bedIntervals mode flatten

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
    }' ${intervals}
    """
  else
    """
    awk -v FS="[:-]" '{
      chunk_fname = sprintf("%s_%d-%d.bed", \$1, \$2, \$3);
      printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > chunk_fname
    }' ${intervals}
    """
}

bedIntervals = bedIntervals
  .map { intervalFile ->
    def duration = 0.0
    for (line in intervalFile.readLines()) {
      final fields = line.split('\t')
      if (fields.size() >= 5) duration += fields[4].toFloat()
      else {
        start = fields[1].toInteger()
        end = fields[2].toInteger()
        duration += (end - start) / params.nucleotidesPerSecond
      }
    }
    [duration, intervalFile]
  }.toSortedList({ a, b -> b[0] <=> a[0] })
  .flatten().collate(2)
  .map{duration, intervalFile -> intervalFile}

if (params.verbose) bedIntervals = bedIntervals.view {
  "  Interv: ${it.baseName}"
}

bamsAll = bamsNormal.join(bamsTumour)

// Manta, Strelka (no splitting by intervals)
(bamsForManta, bamsForStrelka, bamsForStrelkaBP, bamsForPurple, bamsAll) = bamsAll.into(5)

// MuTect, VarDict (partitioning by intervals)
bamsTumourNormalIntervals = bamsAll.combine(bedIntervals)  // cartesian product of 2 channels
(bamsForMutect, bamsForVardict) = bamsTumourNormalIntervals.into(2)

// This will give as a list of unfiltered calls for MuTect.
process RunMutect {
  tag {idSampleTumour + "_vs_" + idSampleNormal + "-" + intervalBed.baseName}

  input:
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal),
                   idSampleTumour, file(bamTumour), file(baiTumour), file(intervalBed) from bamsForMutect
    set file(genomeFasta), file(genomeIndex), file(genomeDict) from Channel.value([genomeFasta, genomeIndex, genomeDict])

  output:
    set val("MuTect"), idPatient, idSampleNormal, idSampleTumour, file("${intervalBed.baseName}_${idSampleTumour}_vs_${idSampleNormal}.vcf") \
        into mutectOutput

  when: !params.onlyQC

  script:
  """ \
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
  Mutect2 \
  -R ${genomeFasta}\
  -I ${bamTumour}  -tumor ${idSampleTumour} \
  -I ${bamNormal} -normal ${idSampleNormal} \
  -L ${intervalBed} \
  -O ${intervalBed.baseName}_${idSampleTumour}_vs_${idSampleNormal}.vcf \
  """
}

mutectOutput = mutectOutput.groupTuple(by:[0,1,2,3])

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
//    set val("FreeBayes"), idPatient, idSampleNormal, idSampleTumour, file("${intervalBed.baseName}_${idSampleTumour}_vs_${idSampleNormal}.vcf") \
//        into freebayesOutput
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
//
//process RunVarDict {
//  tag {idSampleTumour + "_vs_" + idSampleNormal + "-" + intervalBed.baseName}
//
//  input:
//    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumour, file(bamTumour), file(baiTumour), file(intervalBed) \
//        from bamsForVardict
//    file(genomeFasta) from Channel.value(genomeFasta)
//    file(genomeIndex) from Channel.value(genomeIndex)
//
//  output:
//    set val("VarDict"), idPatient, idSampleNormal, idSampleTumour, file("${intervalBed.baseName}_${idSampleTumour}_vs_${idSampleNormal}.vcf") \
//        into vardictOutput
//
//  when: !params.onlyQC
//
//  script:
//  tmpDir = file("tmp").mkdir()
//  """ \
//  unset JAVA_HOME && \
//  export VAR_DICT_OPTS='-Xms750m -Xmx3000m -XX:+UseSerialGC -Djava.io.tmpdir=${tmpDir} && \
//  vardict-java -G ${genomeFasta} \
//  -N ${idSampleNormal} \
//  -b "${bamTumour}|${bamNormal}" \
//  -c 1 -S 2 -E 3 -g 4 --nosv --deldupvar -Q 10 -F 0x700 -f 0.1 \
//  ${intervalBed} \
//  | awk 'NF>=48' \
//  | testsomatic.R \
//  | var2vcf_paired.pl -P 0.9 -m 4.25 -f 0.01 -M  -N "${idSampleTumour}|${idSampleNormal}" \
//  | bcftools filter -m '+' -s 'REJECT' -e 'STATUS !~ ".*Somatic"' \
//  2> /dev/null \
//  | bcftools filter --soft-filter 'LowFreqBias' --mode '+' -e 'FORMAT/AF[0] < 0.02 && FORMAT/VD[0] < 30 && FORMAT/SBF[0] < 0.1 && FORMAT/NM[0] >= 2.0' \
//  | bcftools filter -i 'QUAL >= 0' \
//  | sed 's/\\.*Somatic\\/Somatic/' \
//  | sed 's/REJECT,Description=".*">/REJECT,Description="Not Somatic via VarDict">/' \
//  | awk -F\$'\t' -v OFS='\t' '{if (\$0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", \$4) } {print}' \
//  | awk -F\$'\t' -v OFS='\t' '{if (\$0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", \$5) } {print}' \
//  | awk -F\$'\t' -v OFS='\t' '\$1!~/^#/ && \$4 == \$5 {next} {print}' \
//  | vcfstreamsort \
//  | bgzip -c > ${vardictOutput} \
//  """
//}
//
//
//process RunVarDictGermline {
//  when: false
//
//  script:
//  tmpDir = file("tmp").mkdir()
//  """ \
//  unset JAVA_HOME && \
//  export VAR_DICT_OPTS='-Xms750m -Xmx3000m -XX:+UseSerialGC -Djava.io.tmpdir=${tmpDir} && \
//  vardict-java -G ${genomeFasta} \
//  -N ${idSampleNormal} \
//  -b ${bamNormal} \
//  -c 1 -S 2 -E 3 -g 4 --nosv --deldupvar -Q 10 -F 0x700 -f 0.1 \
//  ${intervalBed} \
//  | teststrandbias.R \
//  | var2vcf_valid.pl -A -N ${idSampleNormal} -E -f 0.1
//  | bcftools filter -i 'QUAL >= 0' \
//  | bcftools filter --soft-filter 'LowFreqBias' --mode '+' -e 'FORMAT/AF[0] < 0.02 && FORMAT/VD[0] < 30 && INFO/SBF < 0.1 && INFO/NM >= 2.0' \
//  | awk -F\$'\t' -v OFS='\t' '{if (\$0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", \$4) } {print}' \
//  | awk -F\$'\t' -v OFS='\t' '{if (\$0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", \$5) } {print}' \
//  | awk -F\$'\t' -v OFS='\t' '\$1!~/^#/ && \$4 == \$5 {next} {print}' \
//  | vcfstreamsort \
//  | bgzip -c > ${vardictOutput} \
//  """
//}
//
//// we are merging the VCFs that are called separatelly for different intervals
//// so we can have a single sorted VCF containing all the calls for a given caller
//
//vcfsToMerge = mutectOutput.mix(vardictOutput)
//if (params.verbose) vcfsToMerge = vcfsToMerge.view {
//  "VCFs To be merged:\n\
//  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
//  Files : ${it[4].fileName}"
//}
//
//process ConcatVCF {
//  tag {variantCaller + "_" + idSampleTumour + "_vs_" + idSampleNormal}
//
//  publishDir "${params.outDir}/VariantCalling/${idPatient}/${"$variantCaller"}", mode: params.publishDirMode
//
//  input:
//    set variantCaller, idPatient, idSampleNormal, idSampleTumour, file(vcFiles) from vcfsToMerge
//    file(genomeIndex) from Channel.value(genomeIndex)
//    file(targetBED) from Channel.value(params.targetBED ? file(params.targetBED) : "null")
//
//  output:
//    // we have this funny *_* pattern to avoid copying the raw calls to publishdir
//    set variantCaller, idPatient, idSampleNormal, idSampleTumour, file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi") into vcfConcatenated
//    // TODO DRY with ConcatVCF
//
//  when: !params.onlyQC
//
//  script:
//  outputFile = "${variantCaller}_${idSampleTumour}_vs_${idSampleNormal}.vcf"
//  options = params.targetBED ? "-t ${targetBED}" : ""
//  """
//  concatenateVCFs.sh -i ${genomeIndex} -c ${task.cpus} -o ${outputFile} ${options}
//  """
//}
//
//if (params.verbose) vcfConcatenated = vcfConcatenated.view {
//  "Variant Calling output:\n\
//  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
//  File  : ${it[4].fileName}"
//}
//
//process RunStrelka {
//  tag {idSampleTumour + "_vs_" + idSampleNormal}
//
//  publishDir "${params.outDir}/VariantCalling/${idPatient}/Strelka", mode: params.publishDirMode
//
//  input:
//    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumour, file(bamTumour), file(baiTumour) \
//      from bamsForStrelka
//    file(targetBED) from Channel.value(params.targetBED ? file(params.targetBED) : "null")
//    set file(genomeFasta), file(genomeIndex), file(genomeDict) from Channel.value([
//      genomeFasta,
//      genomeIndex,
//      genomeDict
//    ])
//
//  output:
//    set val("Strelka"), idPatient, idSampleNormal, idSampleTumour, file("*.vcf.gz"), file("*.vcf.gz.tbi") into strelkaOutput
//
//  when: !params.onlyQC
//
//  script:
//  beforeScript = params.targetBED ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
//  options = params.targetBED ? "--exome --callRegions call_targets.bed.gz" : ""
//  """
//  ${beforeScript}
//  configureStrelkaSomaticWorkflow.py \
//  --Tumour ${bamTumour} \
//  --normal ${bamNormal} \
//  --referenceFasta ${genomeFasta} \
//  ${options} \
//  --runDir Strelka
//
//  python Strelka/runWorkflow.py -m local -j ${task.cpus}
//  mv Strelka/results/variants/somatic.indels.vcf.gz Strelka_${idSampleTumour}_vs_${idSampleNormal}_somatic_indels.vcf.gz
//  mv Strelka/results/variants/somatic.indels.vcf.gz.tbi Strelka_${idSampleTumour}_vs_${idSampleNormal}_somatic_indels.vcf.gz.tbi
//  mv Strelka/results/variants/somatic.snvs.vcf.gz Strelka_${idSampleTumour}_vs_${idSampleNormal}_somatic_snvs.vcf.gz
//  mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi Strelka_${idSampleTumour}_vs_${idSampleNormal}_somatic_snvs.vcf.gz.tbi
//  """
//}
//
//if (params.verbose) strelkaOutput = strelkaOutput.view {
//  "Variant Calling output:\n\
//  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
//  Files : ${it[4].fileName}\n\
//  Index : ${it[5].fileName}"
//}
//
//process RunManta {
//  tag {idSampleTumour + "_vs_" + idSampleNormal}
//
//  publishDir "${params.outDir}/VariantCalling/${idPatient}/Manta", mode: params.publishDirMode
//
//  input:
//    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumour, file(bamTumour), file(baiTumour) \
//      from bamsForManta
//    file(targetBED) from Channel.value(params.targetBED ? file(params.targetBED) : "null")
//    set file(genomeFasta), file(genomeIndex) from Channel.value([
//      genomeFasta,
//      genomeIndex
//    ])
//
//  output:
//    set val("Manta"), idPatient, idSampleNormal, idSampleTumour, file("*.vcf.gz"), file("*.vcf.gz.tbi") into mantaOutput
//    set idPatient, idSampleNormal, idSampleTumour, file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi")
//        into mantaToStrelka
//
//  when: !params.strelkaBP && !params.onlyQC
//
//  script:
//  beforeScript = params.targetBED ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
//  options = params.targetBED ? "--exome --callRegions call_targets.bed.gz" : ""
//  """
//  ${beforeScript}
//  configManta.py \
//  --normalBam ${bamNormal} \
//  --tumorBam ${bamTumour} \
//  --reference ${genomeFasta} \
//  ${options} \
//  --runDir Manta
//
//  python Manta/runWorkflow.py -m local -j ${task.cpus}
//
//  mv Manta/results/variants/candidateSmallIndels.vcf.gz \
//    Manta_${idSampleTumour}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz
//  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
//    Manta_${idSampleTumour}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz.tbi
//  mv Manta/results/variants/candidateSV.vcf.gz \
//    Manta_${idSampleTumour}_vs_${idSampleNormal}.candidateSV.vcf.gz
//  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
//    Manta_${idSampleTumour}_vs_${idSampleNormal}.candidateSV.vcf.gz.tbi
//  mv Manta/results/variants/diploidSV.vcf.gz \
//    Manta_${idSampleTumour}_vs_${idSampleNormal}.diploidSV.vcf.gz
//  mv Manta/results/variants/diploidSV.vcf.gz.tbi \
//    Manta_${idSampleTumour}_vs_${idSampleNormal}.diploidSV.vcf.gz.tbi
//  mv Manta/results/variants/somaticSV.vcf.gz \
//    Manta_${idSampleTumour}_vs_${idSampleNormal}.somaticSV.vcf.gz
//  mv Manta/results/variants/somaticSV.vcf.gz.tbi \
//    Manta_${idSampleTumour}_vs_${idSampleNormal}.somaticSV.vcf.gz.tbi
//  """
//}
//
//if (params.verbose) mantaOutput = mantaOutput.view {
//  "Variant Calling output:\n\
//  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
//  Files : ${it[4].fileName}\n\
//  Index : ${it[5].fileName}"
//}
//
//// Running Strelka Best Practice with Manta indel candidates
//// For easier joining, remaping channels to idPatient, idSampleNormal, idSampleTumour...
//
//bamsForStrelkaBP = bamsForStrelkaBP.map {
//   idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumour, bamTumour, baiTumour ->
//  [idPatientNormal, idSampleNormal, idSampleTumour, bamNormal, baiNormal, bamTumour, baiTumour]
//}.join(mantaToStrelka, by:[0,1,2]).map {
//   idPatientNormal, idSampleNormal, idSampleTumour, bamNormal, baiNormal, bamTumour, baiTumour, mantaCSI, mantaCSIi ->
//  [idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumour, bamTumour, baiTumour, mantaCSI, mantaCSIi]
//}
//
//process RunStrelkaBP {
//  tag {idSampleTumour + "_vs_" + idSampleNormal}
//
//  publishDir "${params.outDir}/VariantCalling/${idPatient}/Strelka", mode: params.publishDirMode
//
//  input:
//    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumour, file(bamTumour), file(baiTumour),
//        file(mantaCSI), file(mantaCSIi) from bamsForStrelkaBP
//    file(targetBED) from Channel.value(params.targetBED ? file(params.targetBED) : "null")
//    set file(genomeFasta), file(genomeIndex), file(genomeDict) from Channel.value([
//      genomeFasta,
//      genomeIndex,
//      genomeDict
//    ])
//
//  output:
//    set val("Strelka"), idPatient, idSampleNormal, idSampleTumour, file("*.vcf.gz"), file("*.vcf.gz.tbi") into strelkaBPOutput
//
//  when: params.strelkaBP && !params.onlyQC
//
//  script:
//  beforeScript = params.targetBED ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
//  options = params.targetBED ? "--exome --callRegions call_targets.bed.gz" : ""
//  """ \
//  ${beforeScript}
//  configureStrelkaSomaticWorkflow.py \
//  --tumor ${bamTumour} \
//  --normal ${bamNormal} \
//  --referenceFasta ${genomeFasta} \
//  --indelCandidates ${mantaCSI} \
//  ${options} \
//  --runDir Strelka
//
//  python Strelka/runWorkflow.py -m local -j ${task.cpus}
//
//  mv Strelka/results/variants/somatic.indels.vcf.gz \
//    StrelkaBP_${idSampleTumour}_vs_${idSampleNormal}_somatic_indels.vcf.gz
//  mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
//    StrelkaBP_${idSampleTumour}_vs_${idSampleNormal}_somatic_indels.vcf.gz.tbi
//  mv Strelka/results/variants/somatic.snvs.vcf.gz \
//    StrelkaBP_${idSampleTumour}_vs_${idSampleNormal}_somatic_snvs.vcf.gz
//  mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
//    StrelkaBP_${idSampleTumour}_vs_${idSampleNormal}_somatic_snvs.vcf.gz.tbi
//  """
//}
//
//if (params.verbose) strelkaBPOutput = strelkaBPOutput.view {
//  "Variant Calling output:\n\
//  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
//  Files : ${it[4].fileName}\n\
//  Index : ${it[5].fileName}"
//}

//// TODO: wrap PURPLE into a single tool. Or better not? Easier to rerun parts
//process RunAmber {
//  tag {idSampleTumour + "_vs_" + idSampleNormal}
//
//  input:
//    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumour, file(bamTumour), file(baiTumour)
//      from bamsForPurple
//    set file(genomeFasta), file(genomeIndex), file(genomeDict), file(purpleHet) from Channel.value([
//      genomeFasta,
//      genomeIndex,
//      genomeDict,
//      purpleHet
//    ])
//
//  output:
//     set idPatient, status, idSample, '{batch}.amber.baf', '{batch}.amber.baf.pcf' into amberOutput
//
//  when: !params.onlyQC
//
//  script:
//  """
//  AMBER
//  -tumor {wildcards.batch}
//  -tumor_bam {input.tumor_bam}
//  -reference {params.normal_name}
//  -reference_bam {input.normal_bam}
//  -ref_genome {input.ref_fa}
//  -bed {input.snp_bed}
//  -threads {threads}
//  -output_dir {params.outdir} 2>&1 | tee {log}
//  """
//}
//
//rule purple_amber:
//    input:
//        tumor_bam  = lambda wc: batch_by_name[wc.batch].tumor.bam,
//        normal_bam = lambda wc: batch_by_name[wc.batch].normal.bam,
//        snp_bed = hpc.get_ref_file(run.genome_build, 'purple_het'),
//        ref_fa = hpc.get_ref_file(run.genome_build, 'fa'),
//    output:
//        'work/{batch}/purple/amber/{batch}.amber.baf',
//        'work/{batch}/purple/amber/{batch}.amber.baf.pcf',
//    params:
//        normal_name = lambda wc: batch_by_name[wc.batch].normal.name,
//        outdir = 'work/{batch}/purple/amber',
//        jar = join(package_path(), 'jars', 'amber-2.3.jar'),
//        xms = 5000,
//        xmx = purple_mem,
//    log:
//        'log/purple/{batch}/{batch}.amber.log',
//    benchmark:
//        'benchmarks/{batch}/purple/{batch}-amber.tsv'
//    resources:
//        mem_mb = purple_mem,
//    threads:
//        threads_per_batch
//    shell:
//        conda_cmd.format('purple') +

//(strelkaIndels, strelkaSNVS) = strelkaOutput.into(2)
//(mantaSomaticSV, mantaDiploidSV) = mantaOutput.into(2)

//vcfForQC = Channel.empty().mix(
//  vcfConcatenated.map {
//    variantcaller, idPatient, idSampleNormal, idSampleTumour, vcf, tbi ->
//    [variantcaller, vcf]
//  },
//  mantaDiploidSV.map {
//    variantcaller, idPatient, idSampleNormal, idSampleTumour, vcf, tbi ->
//    [variantcaller, vcf[2]]
//  },
//  mantaSomaticSV.map {
//    variantcaller, idPatient, idSampleNormal, idSampleTumour, vcf, tbi ->
//    [variantcaller, vcf[3]]
//  },
//  singleMantaOutput.map {
//    variantcaller, idPatient, idSample, vcf, tbi ->
//    [variantcaller, vcf[2]]
//  },
//  strelkaIndels.map {
//    variantcaller, idPatient, idSampleNormal, idSampleTumour, vcf, tbi ->
//    [variantcaller, vcf[0]]
//  },
//  strelkaSNVS.map {
//    variantcaller, idPatient, idSampleNormal, idSampleTumour, vcf, tbi ->
//    [variantcaller, vcf[1]]
//  })
//
//(vcfForBCFtools, vcfForVCFtools) = vcfForQC.into(2)
//
//process RunBcftoolsStats {
//  tag {vcf}
//
//  publishDir "${params.outDir}/Reports/BCFToolsStats", mode: params.publishDirMode
//
//  input:
//    set variantCaller, file(vcf) from vcfForBCFtools
//
//  output:
//    file ("${vcf.simpleName}.bcf.tools.stats.out") into bcfReport
//
//  when: !params.noReports
//
//  script: QC.bcftools(vcf)
//}
//
//if (params.verbose) bcfReport = bcfReport.view {
//  "BCFTools stats report:\n\
//  File  : [${it.fileName}]"
//}
//
//bcfReport.close()
//
//process RunVcftools {
//  tag {vcf}
//
//  publishDir "${params.outDir}/Reports/VCFTools", mode: params.publishDirMode
//
//  input:
//    set variantCaller, file(vcf) from vcfForVCFtools
//
//  output:
//    file ("${vcf.simpleName}.*") into vcfReport
//
//  when: !params.noReports
//
//  script: QC.vcftools(vcf)
//}
//
//if (params.verbose) vcfReport = vcfReport.view {
//  "VCFTools stats report:\n\
//  File  : [${it.fileName}]"
//}
//
//vcfReport.close()
//
//process GetVersionAlleleCount {
//  publishDir "${params.outDir}/Reports/ToolsVersion", mode: params.publishDirMode
//  output: file("v_*.txt")
//  when: 'ascat' in tools && !params.onlyQC
//
//  script:
//  """
//  alleleCounter --version > v_allelecount.txt
//  """
//}
//
//process GetVersionASCAT {
//  publishDir "${params.outDir}/Reports/ToolsVersion", mode: params.publishDirMode
//  output: file("v_*.txt")
//  when: 'ascat' in tools && !params.onlyQC
//
//  script:
//  """
//  R --version > v_r.txt
//  cat ${baseDir}/scripts/ascat.R | grep "ASCAT version" > v_ascat.txt
//  """
//}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def helpMessage() {
  // Display help message
  log.info "UMCCR Cancer Analysis Workflow"
  log.info "    Usage:"
  log.info "       nextflow run variants.nf --sample <file.tsv> --genome <Genome>"
  log.info ""
  log.info "    --sample <file.tsv>"
  log.info "       Specify a TSV file containing paths to sample files."
  log.info "    --test"
  log.info "       Use a test sample."
  log.info "    --noReports"
  log.info "       Disable QC tools and MultiQC to generate a HTML report"
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

def minimalInformationMessage() {
  // Minimal information message
  log.info "Command line: " + workflow.commandLine
  log.info "Profile     : " + workflow.profile
  log.info "Project dir : " + workflow.projectDir
  log.info "Launch dir  : " + workflow.launchDir
  log.info "Work dir    : " + workflow.workDir
  log.info "Executor    : " + config.process.executor
  log.info "Out dir     : " + params.outDir
  log.info "Input path  : " + tsvFile
  log.info "Genome      : " + params.genome
  log.info "Genomes dir : " + params.genome_base
  log.info "Target BED  : " + params.targetBED
  log.info "Containers"
  if (params.containsKey("repository"))
    log.info "  Repository   : " + params.repository
  if (params.containsKey("containerPath"))
    log.info "  ContainerPath: " + params.containerPath
  if (params.containsKey("tag"))
    log.info "  Tag          : " + params.tag
  log.info "Reference files used:"
  log.info "  fasta       :\n\t" + genomeFasta
  log.info "  dbsnp       :\n\t" + dbsnp
  log.info "  intervals   :\n\t" + intervals
}

workflow.onComplete {
  // Display complete message
  this.minimalInformationMessage()
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
  log.info "Error report: " + (workflow.errorReport ?: '-')
}

workflow.onError {
  // Display error message
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
