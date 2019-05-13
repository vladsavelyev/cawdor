#!/usr/bin/env nextflow

if (params.help) exit 0, helpMessage()

genomeFasta = Utils.findRefFile(params, 'genomeFasta')
genomeIndex = Utils.findRefFile(params, 'genomeIndex')
genomeDict  = Utils.findRefFile(params, 'genomeDict')
bwaIndex    = Utils.findRefFile(params, 'bwaIndex')
if (![genomeFasta, genomeIndex, genomeDict, bwaIndex].every())
  exit 1, "Missing reference files for alignment in ${params.genomes_base} for genome ${params.genome}. " +
      "See --help for more information"

if (!params.containsKey("input")) exit 1, 'Please define which samples to work on by providing the --input option'
inputPath = params.input

//inputPath.eachLine { line ->
//  List row = line.split(/\t/)
//  def sname = row[0]
//  def path = row[1]
//  def phenotype = null
//  if (row.size() > 2) {
//    phenotype = row[2]
//  }
//  if (path == "notfound" || !patb) return
//  if (file(path).isFile()) {
//    println "Found file ${path}"
//    bams << [sname, path, null, phenotype]
//  } else {
//    println "Looking at folder ${path}/*_R1*.*"
//    def found = file("${path}/*_R1*.*")
//    assert found.size() == 1, "R1 files in ${path}: ${found}"
//    def file1 = found[0]
//    println "   found ${file1}"
//    def file2 = file(file1.toString().replace('_R1', '_R2'))
//    fastqs << [sname, file1, file2]
//  }
//}
List inputs = Utils.extractSamples(inputPath)
Utils.checkTumorNormal(inputs)
def inputChannel = Channel.from(inputs)

if (!params.targetDepth) params.targetDepth = 4

Utils.startMessage(log, workflow, config, params, inputPath)

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

fastqsChannel = Channel.create()
bamsChannel = Channel.create()
inputChannel.choice(fastqsChannel, bamsChannel) { it.containsKey('fastq1') ? 0 : 1 }

if (params.verbose) fastqsChannel = fastqsChannel.view {
  "Input FastQ: sample: ${it.sample}, file: ${it.fastq1}"
}

if (params.verbose) bamsChannel = bamsChannel.view {
  "Input BAM: sample: ${it.sample}, file: ${it.bam}"
}

fastqsChannel = fastqsChannel.map {
  it -> [it.sample, it.status, it.patient, file(it.fastq1), file(it.fastq2)]
}

bamsChannel = bamsChannel.map {
  it -> [it.sample, it.status, it.patient, file(it.bam), file(it.bai ? it.bai : "null")]
}

//inputChannel = inputChannel
//  .map {
//      idPatient, status, idSample, idLane, file1, file2 ->
//      [idLane, idPatient, status, idSample, file1]
//  }
//  .groupBy {
//    dPatient, status, idSample, idLane, file1, file2 ->
//    idLane
//  }
//  .group(by: 0)
//  .map { idLane, samples -> samples[0] }

//if (params.verbose) inputChannel = inputChannel.view {
//  "Files to downsample: ${it}"
//}

def targDepth = params.targetDepth
def gSize = 3_234_830_000
def readLength = 150
def readsNum = Math.round((gSize / readLength) * targDepth)
if (params.targetDepth > 0) {
  if (params.verbose) {
    println "Downsampling to ${readsNum} reads, which is ${readsNum * 4} lines (unless the file is smaller)"
  }
}

process MapReads {
  tag {sample}

  publishDir "${params.outDir}/Align", mode: 'link'

  input:
  set sample, status, patient, file(file1), file(file2) from fastqsChannel
  set file(genomeFasta), file(bwaIndex) from Channel.value([genomeFasta, bwaIndex])

  output:
  set sample, status, patient, file("*.bam") into mappedBam

  script:
  CN = ""
  if (params.containsKey("sequencing_center")) {
    CN = 'CN:' + params.get("sequencing_center") + '\\t'
  }
  readGroup = "@RG\\tID:${sname}\\t${CN}\\tSM:${sample}\\tPL:illumina"
  bwaMemCmd = "bwa mem -K 100000000 -R \"${readGroup}\" -t ${task.cpus} -M ${genomeFasta}"
  sortMem = Math.min(2, task.memory.toGiga())
  sortCmd = "samtools sort -@ ${task.cpus} -m ${sortMem}G -"
  inp1 = file1
  inp2 = file2
  if (params.targetDepth > 0) {
    inp1 = "<(gunzip -c ${file1} | head -n${readsNum * 2})"
    inp2 = "<(gunzip -c ${file2} | head -n${readsNum * 2})"
  }
  """ \
  ${bwaMemCmd} \
  ${inp1} \
  ${inp2} \
  | ${sortCmd} \
  > ${sample}.bam \
  """
}

process MarkDuplicates {
  tag {sample}

  publishDir params.outDir, mode: 'link', saveAs: { "Dedup/${it}" }

  input:
  set sample, status, patient, file(inputBam) from mappedBam

  output:
  set sample, status, patient, file("${idSample}.dedup.bam"), file("${idSample}.dedup.bai") into readyBams

  when: !params.onlyQC

  script:
  """ \
  gatk --java-options "-Xms500m -Xmx${task.memory.toGiga()}g" \
  MarkDuplicates \
  --MAX_RECORDS_IN_RAM 50000 \
  --INPUT ${inputBam} \
  --TMP_DIR tmp \
  --ASSUME_SORT_ORDER coordinate \
  --CREATE_INDEX true \
  --OUTPUT ${sample}.dedup.bam \
  """
}

readyBams = readyBams.mix(bamsChannel)

if (params.verbose) readyBams = readyBams.view {
  "Ready BAMs: ${it}"
}

if (params.containsKey("snps")) {
  readyBams.into { readyBamsForMutect; readyBams }

  process RunMutect {
    tag {sample}

    publishDir params.outDir, mode: 'link', saveAs: { "Mutect/${it}" }

    input:
    set sample, status, patient, file(bam), file(bai) from readyBamsForMutect
    set file(genomeFasta), file(genomeIndex), file(genomeDict) from Channel.value([genomeFasta, genomeIndex, genomeDict])
    file(snps) from Channel.value(file(params.snps))

    output:
    set sample, status, patient, file("*.vcf.gz"), file("*.vcf.gz.tbi") into mutectOutput

    script:
    outFile = "${sample}.vcf.gz"
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    Mutect2 \
    -R ${genomeFasta}\
    -I ${bam} -tumor ${idSampleTumour} \
    -L ${snps} \
    -O ${outFile} \
    """
  }
}

process RunConpairPileup {
  tag {sample}

  input:
  set sample, status, patient, file(bam), file(bai) from readyBams
  set file(genomeFasta), file(genomeIndex), file(genomeDict) from Channel.value([genomeFasta, genomeIndex, genomeDict])

  output:
  set sample, status, patient, file('*.pileup') into conpairPileup

  script:
  def idxCmd = ""
  if (bai.name  == "null")
    idxCmd = "samtools index ${bam}"
  """
  ${idxCmd}
  run_gatk_pileup_for_sample.py \
  -B ${bam} \
  -O ${sample}.pileup \
  -g ${params.genome} \
  --reference ${genomeFasta}\
  """
}

normals = Channel.create()
tumours = Channel.create()
conpairPileup.choice(normals, tumours) { it[1] }
normals = normals.map { sample, status, patient, pileup -> [patient, sample, pileup] }
tumours = tumours.map { sample, status, patient, pileup -> [patient, sample, pileup] }
normals.join(tumours).into { patientsConcordance; patientsContamination }

process RunConcordance {
  tag {patient}

  publishDir params.outDir

  input:
  set patient, nSample, nPileup, tSample, tPileup from patientsConcordance

  output:
  set patient, file('*.txt') into concordance

  script:
  """
  verify_concordance.py \
   -T ${input.tPileup} \
   -N ${input.nPileupp} \
   -g ${params.genome} \
   -O ${patient}.txt \
   -H \
   -C 0
  """
}

process RunContamination {
  tag {patient}

  publishDir params.outDir

  input:
  set patient, nSample, nPileup, tSample, tPileup from patientsContamination

  output:
  set patient, file('*_normal_cont.txt'), file('*_tumour_cont.txt') into contamination

  script:
  """
  estimate_tumor_normal_contamination.py \
   -T ${input.tPileup} \
   -N ${input.nPileup} \
   -g ${params.genome} \
   | tee \
   >(grep Tumor | sed 's/Tumor sample/Sample/' > ${tSample}_normal_cont.txt \
   | grep Normal | sed 's/Normal sample/Sample/' > ${nSample}_tumour_cont.txt
"""
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
  log.info "       nextflow run contamination.nf --input <file.tsv> --outDir <Dir> --ptfile <path-to-pt-file>"
  log.info ""
  log.info "    --input <file.tsv>"
  log.info "       Specify a TSV file in form of patient\\tphenotype\\tsample\\bam"
  log.info "    --targetDepth N"
  log.info "       depth to downsample the fastq file (default 4)"
  log.info "    --snps <file.bed>"
  log.info "       call at these snps with Mutect2"
  log.info "    --dedup"
  log.info "       deduplicate BAM file"
  log.info "    --help"
  log.info "       show this help"
  log.info "    --verbose"
  log.info "       Adds more verbosity to workflow"
}

workflow.onError {
  // Display error message
  Utils.nextflowMessage(log, workflow)
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}

