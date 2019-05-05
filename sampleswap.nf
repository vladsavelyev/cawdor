#!/usr/bin/env nextflow

if (params.help) exit 0, helpMessage()

if (!Utils.checkExactlyOne([params.containsKey("fastqs")]))
  exit 1, 'Please define which samples to work on by providing the --fastqs option'

def inputPath = file(params.fastqs)
List fastqs = []
inputPath.eachLine { line ->
  List row = line.split(/\t/)
  sname = row[0]
  folder = row[1]
  if (folder == "notfound" || !folder) return
  println "Looking at folder ${folder}/*_R1*.*"
  def found = file("${folder}/*_R1*.*")
  assert found.size() == 1, "R1 files in ${folder}: ${found}"
  println "   found ${found[0]}"
  fastqs << [sname, found[0]]
}
def inputChannel = Channel.from(fastqs)

Utils.startMessage(log, workflow, config, params, inputPath)

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

if (params.verbose) inputChannel = inputChannel.view {
  "Input files: sample: ${it[0]}, fastq: ${it[1]}"
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

def targDepth = 4
def gSize = 3_234_830_000
def readLength = 150
def readsNum = Math.round((gSize / readLength) * targDepth)
if (params.verbose)
  println "Downsampling to ${readsNum} reads, which is ${readsNum * 4} lines (unless the file is smaller)"

process DownsampleFastq {
  tag {sname}

  input:
  set sname, file(file1) from inputChannel

  output:
  set sname, file("*.ds.fq") into downsampledFastq

  script:
  """
  gunzip -c ${file1} | head -n${readsNum * 4} > ${sname}.ds.fq
  """
}

if (params.verbose) downsampledFastq = downsampledFastq.view {
  "Downsampled: ${it}"
}

process RunNGSCheckMateFastq {
  tag {sname}

  publishDir "${params.outDir}", mode: 'link'

  input:
  set sname, file(fastq) from downsampledFastq
  file ptfile from Channel.value(file(params.ptfile))

  output:
  file "*.ncm" into ngsCheckMateFastq

  script:
  assert ptfile.name.endsWith(".pt")
  """
  ngscheckmate_fastq \
  -1 ${fastq} \
  ${ptfile} \
  -p ${task.cpus} \
  > ${sname}.ncm
  """
}

if (params.verbose) ngsCheckMateFastq = ngsCheckMateFastq.view {
  "NCM files: ${it}"
}

def ncmFiles = ngsCheckMateFastq.toList()

process RunNcnVaf {

  publishDir "${params.outDir}", mode: 'link'

  input:
  file(ncmFiles) from ncmFiles

  output:
  set "output_all.txt", "output_corr_matrix.txt", "r_script.r" into ncmResult

  """
  vaf_ncm.py -f -I . -O .
  """
}

if (params.verbose) ncmResult = ncmResult.view {
  "NCM result: ${it[0]}"
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
  log.info "       nextflow run sampleswap.nf --fastqs <file.tsv> --outDir <Dir> --ptfile <path-to-pt-file>"
  log.info ""
  log.info "    --fastqs <file.tsv>"
  log.info "       Specify a TSV file in form of sampleId\\tfolderWithFastqs"
  log.info "    --ptfile <path.pt>"
  log.info "       NGSCheckMate SNPs PT file"
  log.info "    --help"
  log.info "       show this help"
  log.info "    --verbose"
  log.info "       Adds more verbosity to workflow"
}

workflow.onComplete {
  Utils.endMessage(log, workflow, config, params, inputPath)
}

workflow.onError {
  // Display error message
  Utils.nextflowMessage(log, workflow)
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}

