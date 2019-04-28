#!/usr/bin/env nextflow

if (params.help) exit 0, helpMessage()

if (!Utils.checkExactlyOne([params.test, params.containsKey("samples"), params.containsKey("samplesDir")]))
  exit 1, 'Please define which samples to work on by providing exactly one of the --test, --samples or --samplesDir options'

genomeFasta = Utils.findRefFile(params, 'genomeFasta')
genomeIndex = Utils.findRefFile(params, 'genomeIndex')
genomeDict  = Utils.findRefFile(params, 'genomeDict')
bwaIndex    = Utils.findRefFile(params, 'bwaIndex')
if (![genomeFasta, genomeIndex, genomeDict, bwaIndex].every())
  exit 1, "Missing reference files for alignment. See --help for more information"

inputFiles = Channel.empty()
if (params.containsKey("samples")) {
  inputPath = file(params.samples)
  inputFiles = extractSample(inputPath)
} else if (params.containsKey("samplesDir")) {
  inputFiles = extractFastqFromDir(params.samplesDir)
  (inputFiles, fastqTmp) = inputFiles.into(2)
  fastqTmp.toList().subscribe onNext: {
    if (it.size() == 0) {
      exit 1, "No FASTQ files found in --samplesDir directory '${params.samplesDir}'"
    }
  }
  inputPath = params.samplesDir  // used in the reports
} else exit 1, 'No samples were defined with either --samples input.tsv or --samplesDir /input_dir, see --help'

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

process FastQC {
  tag {idPatient + "-" + idRun}

  publishDir "${params.outDir}/Reports/FastQC/${idRun}", mode: params.publishDirMode

  input:
    set idPatient, status, idSample, idRun, file(inputFile1), file(inputFile2) from inputFilesforFastQC

  output:
    file "*_fastqc.{zip,html}" into fastQCreport

  script:
  inputFiles = Utils.isFq(inputFile1) ? "${inputFile1} ${inputFile2}" : "${inputFile1}"
  """
  fastqc -t 2 -q ${inputFiles}
  """
}

if (params.verbose) fastQCreport = fastQCreport.view {
  "FastQC report:\n\
  Files : [${it[0].fileName}, ${it[1].fileName}]"
}

process MapReads {
  tag {idPatient + "-" + idRun}

  input:
    set idPatient, status, idSample, idRun, file(inputFile1), file(inputFile2) from inputFiles
    set file(genomeFasta), file(bwaIndex) from Channel.value([genomeFasta, bwaIndex])

  output:
    set idPatient, status, idSample, idRun, file("${idRun}.bam") into (mappedBam, mappedBamForQC)

  when: !params.onlyQC

  script:
  CN = ""
  if (params.containsKey("sequencing_center")) CN = "CN:${params.sequencing_center}\\t"
  readGroup = "@RG\\tID:${idRun}\\t${CN}PU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
  // adjust mismatch penalty for tumor samples
  extra = status == 1 ? "-B 3" : ""

  bwaMemCmd = "bwa mem -K 100000000 -p -R \"${readGroup}\" ${extra} -t ${task.cpus} -M ${genomeFasta}"
  sortCmd = "samtools sort -@ ${task.cpus} -m 2G -"

  if (Utils.isFq(inputFile1)) {
    """ \
    ${bwaMemCmd} ${inputFile1} ${inputFile2} \
    | ${sortCmd} \
    > ${idRun}.bam
    """
  } else if (SarekUtils.hasExtension(inputFile1, "bam")) {
    """ \
    samtools sort -n -o -l 1 -@ ${task.cpus} -m ${task.memory.toGiga()}G ${inputFile1} \
    | bedtools bamtofastq -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout \
    | ${bwaMemCmd} \
    | ${sortCmd} \
    > ${idRun}.bam
    """
  }
}

// Sort bam whether they are standalone or should be merged
// Borrowed code from https://github.com/guigolab/chip-nf
singleBam = Channel.create()
groupedBam = Channel.create()
mappedBam.groupTuple(by:[0,1,2])
  .choice(singleBam, groupedBam) {it[2].size() > 1 ? 1 : 0}
singleBam = singleBam.map {
  idPatient, status, idSample, idRun, bam ->
  [idPatient, status, idSample, bam]
}

process MergeBams {
  tag {idPatient + "-" + idSample}

  input:
    set idPatient, status, idSample, idRun, file(bam) from groupedBam

  output:
    set idPatient, status, idSample, file("${idSample}.bam") into mergedBam

  when: !params.onlyQC

  script:
  """
  samtools merge --threads ${task.cpus} ${idSample}.bam ${bam}
  """
}

if (params.verbose) singleBam = singleBam.view {
  "Single BAM:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[3].fileName}]"
}

if (params.verbose) mergedBam = mergedBam.view {
  "Merged BAM:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[3].fileName}]"
}

mergedBam = mergedBam.mix(singleBam)

if (params.verbose) mergedBam = mergedBam.view {
  "BAM for MarkDuplicates:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[3].fileName}]"
}

if (params.verbose) mappedBam = mappedBam.view {
  "Mapped BAM (single or to be merged):\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\tRun   : ${it[3]}\n\
  File  : [${it[4].fileName}]"
}

process MarkDuplicates {
  tag {idPatient + "-" + idSample}

  publishDir params.outDir, mode: 'link',
    saveAs: {
      if (it == "${idSample}.dedup_metrics.txt") "Reports/MarkDuplicates/${it}"
      else "Align/${it}"
    }

  input:
    set idPatient, status, idSample, file("${idSample}.bam") from mergedBam

  output:
    set idPatient, idSample, file("${idSample}.dedup.bam"), file("${idSample}.dedup.bai") into duplicateMarkedBams
    set idPatient, status, idSample, val("${idSample}.dedup.bam"), val("${idSample}.dedup.bai") into samplesTsv
    file ("${idSample}.dedup_metrics.txt") into markDuplicatesReport

  when: !params.onlyQC

  script:
  """ \
  gatk --java-options "-Xms500m -Xmx${task.memory.toGiga()}g" \
  MarkDuplicates \
  --MAX_RECORDS_IN_RAM 50000 \
  --INPUT ${idSample}.bam \
  --METRICS_FILE ${idSample}.dedup_metrics.txt \
  --TMP_DIR . \
  --ASSUME_SORT_ORDER coordinate \
  --CREATE_INDEX true \
  --OUTPUT ${idSample}.dedup.bam \
  """
}

// Creating a TSV file to restart from this step
samplesTsv.map { idPatient, status, idSample, bam, bai ->
  "${idPatient}\t${status}\t${idSample}\t${params.outDir}/Align/${bam}\t${params.outDir}/Align/${bai}\n"
}.collectFile(
  name: 'samples.tsv', sort: true, storeDir: "${params.outDir}"
)

(bamForQualimap, bamForSamToolsStats) = duplicateMarkedBams.into(2)

process RunQualimap {
  tag {idPatient + "-" + idSample}

  publishDir "${params.outDir}/Reports/Qualimap", mode: params.publishDirMode

  input:
  set idPatient, idSample, file(bam), file(bai) from bamForQualimap

  output:
  file("${bam.baseName}") into qualimapReport

  when: !params.noQualimap

  script:
  """
  qualimap --java-mem-size=${task.memory.toGiga()}G \
  bamqc \
  -bam ${bam} \
  --paint-chromosome-limits \
  --genome-gc-distr HUMAN \
  -nt ${task.cpus} \
  -skip-duplicated \
  --skip-dup-mode 0 \
  -outdir ${bam.baseName} \
  -outformat HTML
  """
}

if (params.verbose) qualimapReport = qualimapReport.view {
  "Qualimap BamQC report:\n\
  Dir   : [${it.fileName}]"
}

process RunSamtoolsStats {
  tag {idPatient + "-" + idSample}

  publishDir "${params.outDir}/Reports/SamToolsStats", mode: params.publishDirMode

  input:
    set idPatient, idSample, file(bam), file(bai) from bamForSamToolsStats

  output:
    file ("${bam}.stats.txt") into samtoolsStatsReport

  script:
  """ \
  samtools stats ${bam} > ${bam}.stats.txt \
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
    else if (file1.isDirectory()) {
      row = extractFastqFromDir(file1)
      file1 = row[4]
      file2 = row.size() > 5 ? row[5] : file("null")
    } else {
      "No recognisable extention for input file: ${file1}"
    }

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
      for (path1 in file("${sampleDir}/**_R1_*.fastq.gz")) {
        assert path1.getName().contains('_R1_')
        path2 = file(path1.toString().replace('_R1_', '_R2_'))
        if (!path2.exists()) error "Path '${path2}' not found"
        // the last name of the sampleDir is assumed to be a unique sample id
        sampleId = path1.getParent().getName().toString()
        patient = sampleId
        (flowcell, lane) = flowcellLaneFromFastq(path1)
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
  log.info "Command line: " + workflow.commandLine
  log.info "Profile     : " + workflow.profile
  log.info "Project dir : " + workflow.projectDir
  log.info "Launch dir  : " + workflow.launchDir
  log.info "Work dir    : " + workflow.workDir
  log.info "Cont Engine : " + workflow.containerEngine
  log.info "Executor    : " + config.process.executor
  log.info "Out dir     : " + params.outDir
  log.info "Input path  : " + inputPath
  log.info "Genome      : " + params.genome
  log.info "Genomes dir : " + params.genome_base
  log.info "Containers"
  if (params.containsKey("repository"))
    log.info "  Repository   : " + params.repository
  if (params.containsKey("containerPath"))
    log.info "  ContainerPath: " + params.containerPath
  if (params.containsKey("tag"))
    log.info "  Tag          : " + params.tag
  log.info "Reference files used:"
  log.info "  fasta       :\n\t" + genomeFasta
  log.info "  bwa indexes :\n\t" + bwaIndex.join(',\n\t')
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

