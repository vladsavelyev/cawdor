import static nextflow.Nextflow.file
import static nextflow.Nextflow.exit
import static nextflow.Nextflow.error
import nextflow.Channel

class Utils {

  static def findRefFile(params, key) {
    def genomes = params.genomes
    def genome = params.genome

    if (!(genome in genomes)) {
      exit 1, "Genome ${genome} not found in configuration. Available: ${genomes}"
    }

    def path = genomes[genome][key]

    if (!path) {
      path = genomes["default"][key]
      if (!path) exit 1, "Missing key ${key} for genome ${genome}"
    }

    path = file(params.genomes_base + "/" + path)
    if (!Utils.checkRefExistence(key, path)) {
      System.err.println "Not found reference file ${key} ${path}"
      path = null
    }
    return path
  }

  // Loop through all the references files to check their existence
  static boolean checkRefExistence(key, fileToCheck) {
    if (fileToCheck instanceof List) return fileToCheck.every{ Utils.checkRefExistence(key, it) }
    def f = file(fileToCheck)
    // this is an expanded wildcard: we can assume all files exist
    if (f instanceof List && f.size() > 0) return true
    else if (!f.exists()) {
      System.err.println "Missing reference files: ${key} ${fileToCheck}"
      return false
    }
    return true
  }

  static boolean checkExactlyOne(list) {
    def n = 0
    list.each{n += it ? 1 : 0}
    return n == 1
  }

  // Return status [0, 1]
  // 0 == Normal, 1 == Tumor
  static def returnStatus(it) {
    if (it == 'T') it = '1'
    if (it == 'N') it = '0'
    if (!(it in ['0', '1'])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it.toInteger()
  }

  // Return file if it exists
  static def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
  }

  // Check if a row has the expected number of item
  static boolean checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
  }

  // Check file extension
  static def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
  }

  // Check if the file extention corresponds to fastq
  static boolean isFq(it) {
    hasExtension(it, 'fastq.gz') || hasExtension(it, 'fq.gz')
  }

  // Check parameter existence
  static boolean checkParameterExistence(it, list) {
    if (!list.contains(it)) {
      println "Unknown parameter: ${it}"
      return false
    }
    return true
  }

  // Compare each parameter with a list of parameters
  static boolean checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
  }

  // Channeling the TSV file containing BAM.
  // Format is: "subject status sample bam bai"
  static def extractBams(tsvFile) {

    Map tumoursByPatient = [:]
    Map normalsByPatient = [:]

    tsvFile.eachLine { line ->
      if (line ==~ /^#.*/) return
      List row = line.split(/\t/)
      if (row.size() == 0) return
      println("TSV row: ${row}")

      if (row.size() < 4) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
      def idPatient = row[0]
      def status    = Utils.returnStatus(row[1])
      def idSample  = row[2]
      def bamFile   = Utils.returnFile(row[3])
      def baiFile
      if (!Utils.hasExtension(bamFile, ".bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
      if (row.size() >= 5) {
        baiFile   = Utils.returnFile(row[4])
        if (!Utils.hasExtension(baiFile, ".bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"
      } else {
        baiFile = bamFile + ".bai"
      }
      if (status == 0) normalsByPatient[idPatient] = [idSample, bamFile, baiFile]
      else             tumoursByPatient[idPatient] = [idSample, bamFile, baiFile]
    }

    Boolean isErr = false
    if (tumoursByPatient.keySet() - normalsByPatient.keySet()) {
      System.err.println "Found tumours without normals with patientIDs: " + \
          "${tumoursByPatient.keySet() - normalsByPatient.keySet()}. Check your TSV"
      isErr = true
    }
    if (normalsByPatient.keySet() - tumoursByPatient.keySet()) {
      System.err.println "Found normals without tumours with patientIDs: " + \
          "${normalsByPatient.keySet() - tumoursByPatient.keySet()}. Check your TSV"
      isErr = true
    }
    if (isErr) exit 1

    tumoursByPatient.collect { patId, tumours ->
      [patId] + normalsByPatient[patId] + tumours
    }
  }

  static List extractSamplesFromTSV(tsvFile) {
    // Channeling the TSV file containing FASTQ or BAM
    // Format is: "subject status sample lane fastq1 fastq2"
    // or: "subject status sample lane bam"
    List samplesFromTSV = []

    tsvFile.eachLine { line ->
      if (line ==~ /^#.*/) return
      List row = line.split(/\t/)
      if (row.size() == 0) return
      println("TSV row: ${row}")

      String idPatient = row[0]
      int status       = Utils.returnStatus(row[1])
      String idSample  = row[2]
      String idLane    = "${idSample}.${row[3]}"
      def file1        = Utils.returnFile(row[4])
      def file2        = file("null")

      if (Utils.isFq(file1)) {
        if (row.size() > 5) {
          println("Paths are FastQ files: ${file1} and ${file2}")
          file2 = Utils.returnFile(row[5])
          if (!Utils.isFq(file2)) exit 1, "R2 FastQ file ${file2} has the wrong extension. " +
              "See --help for more information"
        }
        else {
          println("Path ${file1} is a R1 fastq file. Found adjacent R2 file: ${file2}.")
          if (!file1.getName().contains('_R1')) exit 1, "Can't find R2 match for ${file1} as the file name " +
              "doesn't contain _R1. Workaround is to specify both R1 and R2 FastQ files in the TSV. " +
              "See --help for more information"
          file2 = file(file1.toString().replace('_R1', '_R2'))
          if (!file2.exists()) error "Not found R2 FastQ file '${file2}'"
        }
        samplesFromTSV << [idPatient, status, idSample, idLane, file1, file2]
      }
      else if (file1.toString().toLowerCase().endsWith(".bam")) {
        println("Path ${file1} is a BAM file")
        Utils.checkNumberOfItem(row, 5)
        samplesFromTSV << [idPatient, status, idSample, idLane, file1, file2]
      }
      else if (file1.isDirectory()) {
        println("Path ${file1} is a directory, searching it for FastQ files")
        List fastqsFromDir = Utils.extractFastqFromDir(file1)
        fastqsFromDir.each { r ->
          def f1 = Utils.returnFile(r[4])
          def f2 = r.size() > 5 ? Utils.returnFile(r[5]) : file("null")
          samplesFromTSV << [idPatient, status, idSample, r[3], f1, f2]
        }
      } else {
        "No recognisable extention for input file: ${file1}"
      }
    }

    samplesFromTSV
  }

  static List extractFastqFromDir(pattern) {
    // create a channel of FASTQs from a directory pattern such as
    // "my_samples/*/". All samples are considered 'normal'.
    // All FASTQ files in subdirectories are collected and emitted;
    // they must have _R1 and _R2 in their names.

    List res = []
    List dirs = [file(pattern, type: "dir")].flatten()
    dirs.each {
      for (path1 in file("${it}/**_R1*.fastq.gz")) {
        assert path1.getName().contains('_R1')
        def path2 = file(path1.toString().replace('_R1', '_R2'))
        if (!path2.exists()) error "Path '${path2}' not found"
        // the last name of the sampleDir is assumed to be a unique sample id
        String sampleId = path1.getParent().getName().toString()
        String patient = sampleId
        def (String flowcell, int lane) = Utils.flowcellLaneFromFastq(path1)
        int status = 0  // normal (not tumor)
        GString rgId = "${flowcell}.${sampleId}.${lane}"
        res << [patient, status, sampleId, rgId, path1, path2]
      }
    }
    res
  }

  static def flowcellLaneFromFastq(path) {
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

  static def cawdorMsg(log) {
    log.info "Cawdor - UMCCR cancer analysis workflow for DNA and RNA WGS sequencing data"
  }

  static def nextflowMessage(log, workflow) {
    // Nextflow message (version + build)
    log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
  }

  static def startMessage(log, workflow, config, params) {
    // Minimal information message
    log.info "Command line: " + workflow.commandLine
    log.info "Profile     : " + workflow.profile
    log.info "Project dir : " + workflow.projectDir
    log.info "Launch dir  : " + workflow.launchDir
    log.info "Work dir    : " + workflow.workDir
    log.info "Cont engine : " + workflow.containerEngine
    log.info "Executor    : " + config.process.executor
    log.info "Out dir     : " + params.outDir
    log.info "Genome      : " + params.genome
    log.info "Genomes dir : " + params.genomes_base
  }

  static def endMessage(log, workflow, config, params) {
    startMessage(log, workflow, config, params)
    log.info "Completed at: " + workflow.complete
    log.info "Duration    : " + workflow.duration
    log.info "Success     : " + workflow.success
    log.info "Exit status : " + workflow.exitStatus
    log.info "Error report: " + (workflow.errorReport ?: '-')
  }
}