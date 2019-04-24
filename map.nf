#!/usr/bin/env nextflow

import Utils

if (params.help) exit 0, helpMessage()

referenceMap = Utils.defineReferenceMap(params)

if (!Utils.checkExactlyOne([params.test, params.sample, params.sampleDir]))
  exit 1, 'Please define which samples to work on by providing exactly one of the --test, --sample or --sampleDir options'
if (!Utils.checkReferenceMap(referenceMap)) exit 1, 'Missing Reference file(s), see --help for more information'

tsvPath = ''
if (params.sample) tsvPath = params.sample

inputFiles = Channel.empty()
if (tsvPath) {
  tsvFile = file(tsvPath)
  inputFiles = extractSample(tsvFile)
} else if (params.sampleDir) {
  inputFiles = extractFastqFromDir(params.sampleDir)
  (inputFiles, fastqTmp) = inputFiles.into(2)
  fastqTmp.toList().subscribe onNext: {
    if (it.size() == 0) {
      exit 1, "No FASTQ files found in --sampleDir directory '${params.sampleDir}'"
    }
  }
  tsvFile = params.sampleDir  // used in the reports
} else exit 1, 'No sample were defined, see --help'

minimalInformationMessage()

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

// fork inputFiles channel into 2 copies
(inputFiles, inputFilesforFastQC) = inputFiles.into(2)

if (params.verbose) inputFiles = inputFiles.view {
  "Input files:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\tRun   : ${it[3]}\n\
  Files : [${it[4].fileName}, ${it[5].fileName}]"
}

process MapReads {
  tag {idPatient + "-" + idRun}

  input:
    set idPatient, status, idSample, idRun, file(inputFile1), file(inputFile2) from inputFiles
    set file(genomeFasta), file(bwaIndex) from Channel.value([referenceMap.genomeFasta, referenceMap.bwaIndex])

  output:
    set idPatient, status, idSample, idRun, file("${idRun}.bam"), file("${idRun}.bam.bai") into (mappedBam, mappedBamForQC)

  when: !params.onlyQC

  script:
  CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
  readGroup = "@RG\\tID:${idRun}\\t${CN}PU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
  // adjust mismatch penalty for tumor samples
  extra = status == 1 ? "-B 3" : ""

  bwaMemCmd = "bwa mem -K 100000000 -p -R \"${readGroup}\" ${extra} -t ${task.cpus} -M ${genomeFasta}"
  sortCmd = "samtools sort -n -@ ${task.cpus} -m 2G -O sam -"
  dedupCmd = "bamsormadup inputformat=sam threads=${task.cpus} SO=coordinate"

  if (Utils.isFq(inputFile1)) {
    """ \
    ${bwaMemCmd} ${inputFile1} ${inputFile2} \
    | ${sortCmd} \
    | ${dedupCmd} \
    > ${idRun}.bam \
    && samtools index ${idRun}.bam
    """
  } else if (SarekUtils.hasExtension(inputFile1, "bam")) {
    """ \
    samtools sort -n -o -l 1 -@ ${task.cpus} -m ${task.memory.toGiga()}G ${inputFile1} \
    | bedtools bamtofastq -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout \
    | ${bwaMemCmd} \
    | ${sortCmd} \
    | ${dedupCmd} \
    > ${idRun}.bam \
    && samtools index ${idRun}.bam
    """
  }
}

if (params.verbose) mappedBam = mappedBam.view {
  "Mapped BAM (single or to be merged):\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\tRun   : ${it[3]}\n\
  File  : [${it[4].fileName}]"
}

process RunSamtoolsStats {
  tag {idPatient + "-" + idSample}

  publishDir "${params.outDir}/QC/SamToolsStats", mode: params.publishDirMode

  input:
    set idPatient, status, idSample, idRun, file(bam), file(bai) from mappedBamForQC

  output:
    file ("${bam}.stats.txt") into samtoolsStatsReport

  script:
  """
  samtools stats ${bam} > ${bam}.stats.txt
  """
}

if (params.verbose) samtoolsStatsReport = samtoolsStatsReport.view {
  "SAMTools stats report:\n\
  File  : [${it.fileName}]"
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def extractSample(tsvFile) {
  // Channeling the TSV file containing FASTQ or BAM
  // Format is: "subject status sample lane fastq1 fastq2"
  // or: "subject status sample lane bam"
  Channel.from(tsvFile)
  .splitCsv(sep: '\t')
  .map { row ->
    def idPatient  = row[0]
    def status     = Utils.returnStatus(row[1].toInteger())
    def idSample   = row[2]
    def idRun      = row[3]
    def file1      = Utils.returnFile(row[4])
    def file2      = file("null")
    if (Utils.isFq(file1)) {
      Utils.checkNumberOfItem(row, 6)
      file2 = Utils.returnFile(row[5])
      if (!Utils.isFq(file2)) exit 1, "File: ${file2} has the wrong extension. See --help for more information"
    }
    else if (file1.toString().toLowerCase().endsWith(".bam")) {
      Utils.checkNumberOfItem(row, 5)
      if (!Utils.hasExtension(file1, "bam")) exit 1, "File: ${file1} has the wrong extension. See --help for more information"
    }
    else "No recognisable extention for input file: ${file1}"

    [idPatient, status, idSample, idRun, file1, file2]
  }
}

def extractFastqFromDir(pattern) {
  // create a channel of FASTQs from a directory pattern such as
  // "my_samples/*/". All samples are considered 'normal'.
  // All FASTQ files in subdirectories are collected and emitted;
  // they must have _R1_ and _R2_ in their names.

  def fastq = Channel.create()

  // a temporary channel does all the work
  Channel
    .fromPath(pattern, type: 'dir')
    .ifEmpty { error "No directories found matching pattern '${pattern}'" }
    .subscribe onNext: { sampleDir ->
      // the last name of the sampleDir is assumed to be a unique sample id
      sampleId = sampleDir.getFileName().toString()

      for (path1 in file("${sampleDir}/**_R1_*.{fastq,fq}.gz")) {
        assert path1.getName().contains('_R1_')
        path2 = file(path1.toString().replace('_R1_', '_R2_'))
        if (!path2.exists()) error "Path '${path2}' not found"
        (flowcell, lane) = flowcellLaneFromFastq(path1)
        patient = sampleId
        status = 0  // normal (not tumor)
        rgId = "${flowcell}.${sampleId}.${lane}"
        result = [patient, status, sampleId, rgId, path1, path2]
        fastq.bind(result)
      }
  }, onComplete: { fastq.close() }

  fastq
}

def flowcellLaneFromFastq(path) {
  // parse first line of a FASTQ file (optionally gzip-compressed)
  // and return the flowcell id and lane number.
  // expected format:
  // xx:yy:FLOWCELLID:LANE:... (seven fields)
  // or
  // FLOWCELLID:LANE:xx:... (five fields)
  InputStream fileStream = new FileInputStream(path.toFile())
  InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
  Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
  BufferedReader buffered = new BufferedReader(decoder)
  def line = buffered.readLine()
  assert line.startsWith('@')
  line = line.substring(1)
  def fields = line.split(' ')[0].split(':')
  String fcid
  int lane
  if (fields.size() == 7) {
    // CASAVA 1.8+ format
    fcid = fields[2]
    lane = fields[3].toInteger()
  }
  else if (fields.size() == 5) {
    fcid = fields[0]
    lane = fields[1].toInteger()
  }
  [fcid, lane]
}

def helpMessage() {
  // Display help message
  log.info "UMCCR Cancer Analysis Workflow"
  log.info "    Usage:"
  log.info "       nextflow run main.nf --sample <file.tsv> --genome <Genome>"
  log.info "       nextflow run main.nf --sampleDir <Directory> --genome <Genome>"
  log.info ""
  log.info "    --sample <file.tsv>"
  log.info "       Specify a TSV file containing paths to sample files."
  log.info "    --sampleDir <Directoy>"
  log.info "       Specify a directory containing sample files."
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
  log.info "       show this help"
  log.info "    --verbose"
  log.info "       Adds more verbosity to workflow"
}

def minimalInformationMessage() {
  // Minimal information message
  log.info "Command Line: " + workflow.commandLine
  log.info "Profile     : " + workflow.profile
  log.info "Project Dir : " + workflow.projectDir
  log.info "Launch Dir  : " + workflow.launchDir
  log.info "Work Dir    : " + workflow.workDir
  log.info "Cont Engine : " + workflow.containerEngine
  log.info "Out Dir     : " + params.outDir
  log.info "TSV file    : ${tsvFile}"
  log.info "Genome      : " + params.genome
  log.info "Genome_base : " + params.genome_base
  log.info "Containers"
  if (params.repository != "") log.info "  Repository   : " + params.repository
  if (params.containerPath != "") log.info "  ContainerPath: " + params.containerPath
  log.info "  Tag          : " + params.tag
  log.info "Reference files used:"
  log.info "  genome      :\n\t" + referenceMap.genomeFasta
  log.info "\t" + referenceMap.genomeDict
  log.info "\t" + referenceMap.genomeIndex
  log.info "  bwa indexes :\n\t" + referenceMap.bwaIndex.join(',\n\t')
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

