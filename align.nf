#!/usr/bin/env nextflow

if (params.help) exit 0, helpMessage()

if (!Utils.checkExactlyOne([params.test, params.containsKey("samples"), params.containsKey("samplesDir")]))
  exit 1, 'Please define which samples to work on by providing exactly one of the --test, --samples or --samplesDir options'

genomeFasta = Utils.findRefFile(params, 'genomeFasta')
genomeIndex = Utils.findRefFile(params, 'genomeIndex')
genomeDict  = Utils.findRefFile(params, 'genomeDict')
bwaIndex    = Utils.findRefFile(params, 'bwaIndex')
if (![genomeFasta, genomeIndex, genomeDict, bwaIndex].every())
  exit 1, "Missing reference files for alignment in ${params.genomes_base} for genome ${params.genome}. " +
          "See --help for more information"

def inputChannel = Channel.empty()
def inputPath = file("null")
if (params.containsKey("samples")) {
  inputPath = file(params.samples)
  List inputFiles = Utils.extractSamplesFromTSV(inputPath)
  inputChannel = Channel.from(inputFiles)
  (inputChannel, tmp) = inputChannel.into(2)
  tmp.toList().subscribe onNext: {
    if (it.size() == 0) {
      exit 1, "No FASTQ files found in TSV file '${params.samples}'"
    }
  }
} else if (params.containsKey("samplesDir")) {
  inputPath = params.samplesDir  // used in the reports
  List inputFiles = Utils.extractFastqFromDir(params.samplesDir)
  inputChannel = Channel.from(inputFiles)
  (inputChannel, tmp) = inputChannel.into(2)
  tmp.toList().subscribe onNext: {
    if (it.size() == 0) {
      exit 1, "No FASTQ files found in --samplesDir directory '${params.samplesDir}'"
    }
  }
} else {
  exit 1, 'No samples were defined with either --samples input.tsv or --samplesDir /input_dir, see --help'
}

Utils.startMessage(log, workflow, config, params, inputPath)

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

// fork inputFiles channel into 2 copies
(inputChannel, inputFilesforFastQC) = inputChannel.into(2)

if (params.verbose) inputChannel = inputChannel.view {
  "Input files:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\tRun   : ${it[3]}\n\
  Files : [${it[4].fileName}, ${it[5].fileName}]"
}

process RunFastQC {
  tag {idPatient + "-" + idLane}

  publishDir "${params.outDir}/Reports/FastQC/${idLane}", mode: params.publishDirMode

  input:
  set idPatient, status, idSample, idLane, file(inputFile1), file(inputFile2) from inputFilesforFastQC

  output:
  file "*_fastqc.{zip,html}" into fastQCreport

  script:
  def inputs = Utils.isFq(inputFile1) ? "${inputFile1} ${inputFile2}" : "${inputFile1}"
  """
  fastqc -t 2 -q ${inputs}
  """
}

if (params.verbose) fastQCreport = fastQCreport.view {
  "FastQC report:\n\
  Files : [${it[0].fileName}, ${it[1].fileName}]"
}

process MapReads {
  tag {idPatient + "-" + idLane}

  input:
  set idPatient, status, idSample, idLane, file(inputFile1), file(inputFile2) from inputChannel
  set file(genomeFasta), file(bwaIndex) from Channel.value([genomeFasta, bwaIndex])

  output:
  set idPatient, status, idSample, idLane, file("${idLane}.bam") into mappedBam

  when: !params.onlyQC

  script:
  CN = ""
  if (params.containsKey("sequencing_center")) {
    CN = 'CN:' + params.get("sequencing_center") + '\\t'
  }
  readGroup = "@RG\\tID:${idLane}\\t${CN}PU:${idLane}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
  // adjust mismatch penalty for tumor samples
  extra = status == 1 ? "-B 3" : ""
//  bwa mem -K 100000000 -p -R "@RG\tID:D0EN0ACXX111207.normal.4\tPU:D0EN0ACXX111207.normal.4\tSM:normal\tLB:normal\tPL:illumina"  -t 28 -M human_g1k_v37_decoy.small.fasta tiny_n_L004_R1_xxx.fastq.gz tiny_n_L004_R2_xxx.fastq.gz
  bwaMemCmd = "bwa mem -K 100000000 -R \"${readGroup}\" ${extra} -t ${task.cpus} -M ${genomeFasta}"
  sortMem = Math.min(2, task.memory.toGiga())
  sortCmd = "samtools sort -@ ${task.cpus} -m ${sortMem}G -"

  if (Utils.isFq(inputFile1)) {
    """ \
    ${bwaMemCmd} ${inputFile1} ${inputFile2} \
    | ${sortCmd} \
    > ${idLane}.bam \
    """
  } else if (Utils.hasExtension(inputFile1, "bam")) {
    """ \
    samtools sort -n -o -l 1 -@ ${task.cpus} -m ${task.memory.toGiga()}G ${inputFile1} \
    | bedtools bamtofastq -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout \
    | ${bwaMemCmd} -p \
    | ${sortCmd} \
    > ${idLane}.bam \
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
  idPatient, status, idSample, idLane, bam ->
  [idPatient, status, idSample, bam]
}

process MergeBams {
  tag {idPatient + "-" + idSample}

  input:
  set idPatient, status, idSample, idLane, file(bam) from groupedBam

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
  "All per-sample BAMs to be deduplicated:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[3].fileName}]"
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
  --TMP_DIR tmp \
  --ASSUME_SORT_ORDER coordinate \
  --CREATE_INDEX true \
  --OUTPUT ${idSample}.dedup.bam \
  """
}

// Creating a TSV file to restart from this step
samplesTsv.map { idPatient, status, idSample, bam, bai ->
  "${idPatient}\t${status}\t${idSample}\t${params.outDir}/Align/${bam}\t${params.outDir}/Align/${bai}\n"
}.collectFile(
  name: 'samples.tsv', sort: true, storeDir: "${params.outDir}/Align"
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

def helpMessage() {
  // Display help message
  Utils.cawdorMsg(log)
  log.info "    Usage:"
  log.info "       nextflow run align.nf --samples <file.tsv> --outDir <Dir> --genome <Genome>"
  log.info "       nextflow run align.nf --samplesDir <Directory> --outDIr <Dir> --genome <Genome>"
  log.info ""
  log.info "    --samples <file.tsv>"
  log.info "       Specify a TSV file containing paths to sample files."
  log.info "    --samplesDir <Directoy>"
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

workflow.onComplete {
  Utils.endMessage(log, workflow, config, params, inputPath)
}

workflow.onError {
  // Display error message
  Utils.nextflowMessage(log, workflow)
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}

