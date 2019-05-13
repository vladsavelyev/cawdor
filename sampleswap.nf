#!/usr/bin/env nextflow

if (params.help) exit 0, helpMessage()

if (!Utils.checkExactlyOne([params.containsKey("fastqs")]))
  exit 1, 'Please define which samples to work on by providing the --fastqs option'

def inputPath = file(params.fastqs)
List fastqs = []
inputPath.eachLine { line ->
  List row = line.split(/\t/)
  sname = row[0]
  path = row[1]
  if (path == "notfound" || !path) return
  if (file(path).isFile()) {
    println "Found file ${path}"
    fastqs << [sname, file1]
  } else if (file(path).isDirectory()) {
    println "Looking at folder ${path}"
    def found = file("${path}/*_R1*.*")
    if (found.size() > 0)
    assert found.size() == 1, "R1 files in ${path}: ${found}"
    def file1 = found[0]
    println "   found ${file1}"
    fastqs << [sname, file1]
  }
}
def inputChannel = Channel.from(fastqs)

if (!params.targetDepth) params.targetDepth = 4

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

if (params.targetDepth > 0) {
  def targDepth = params.targetDepth
  def gSize = 3_234_830_000
  def readLength = 150
  def readsNum = Math.round((gSize / readLength) * targDepth)
  if (params.verbose) {
    println "Downsampling to ${readsNum} reads, which is ${readsNum * 4} lines (unless the file is smaller)"
  }

  process DownsampleFastq {
    tag {sname}

    validExitStatus 0,141

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

  inputChannel = downsampledFastq
}

process RunNcmFastq {
  tag {sname}

  publishDir "${params.outDir}", mode: 'link'

  input:
  set sname, file(file1) from inputChannel
  file ptfile from Channel.value(file(params.ptfile))

  output:
  file "*.ncm" into ncmFastq

  script:
  assert ptfile.name.endsWith(".pt")
  """
  ngscheckmate_fastq \
  -1 ${file1} \
  ${ptfile} \
  -p ${task.cpus} \
  > ${sname}.ncm
  """
}

if (params.verbose) ncmFastq = ncmFastq.view {
  "NCM files: ${it}"
}

def ncmFiles = ncmFastq.toList()

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
  log.info "    --targetDepth N"
  log.info "       depth to downsample the fastq file (default 4. Set to -1 to skip)"
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

