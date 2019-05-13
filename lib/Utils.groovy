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
    hasExtension(it, '.fastq.gz') || hasExtension(it, '.fq.gz') || hasExtension(it, '.fastq') || hasExtension(it, '.fq')
  }

  static boolean checkBai(bai, bam) {
    // checks if the file has the .bai extension and the matching BAM has the same base name
    Utils.hasExtension(bai, '.bai') && (bai.startsWith(bam -~/\.bam/) || bai.startsWith(bam -~/\.BAM/))
  }

  // Channeling the TSV file containing BAM.
  // Format is: "subject status sample bam bai"
  static List extractSamples(tsvFile) {
    List header = null
    List samples = []

    file(tsvFile).eachLine { line ->
      line = line.trim()
      if (!line) return
      List row = line.split(/\t/)
      println("TSV row: ${row}")
      if (line ==~ /^#.*/) {
        row = (line -~ /^#/).split(/\t/)
        if (row.size() > 0) header = row
      } else {
        Map data = [:]

        if (header) {
          if (row.size() < header.size()) exit 1, "Malformed row in TSV file: ${row}, expected format: ${header}"
        }

        if (header == null) {
          if (row.size() == 1) { header = ['file1'] }
          if (row.size() == 2) {
            header = ['sample', 'file1']
            if (hasExtension(row[1], '.bai') || Utils.isFq(row[0]) && Utils.isFq(row[1]))
              header = ['file1', 'file2']
          }
          if (row.size() == 3) { header = ['sample', 'file1', 'file2'] }
          if (row.size() == 4) { header = ['patient', 'status', 'sample', 'file1'] }
          if (row.size() == 5) { header = ['patient', 'status', 'sample', 'file1', 'file2'] }
        }
        header.eachWithIndex { h, i -> data[h] = row[i] }

        def availFileKeys = ['path', 'file', 'file1', 'fastq', 'fastq1', 'bam', 'dir', 'folder'].toSet()
        def foundFileKeys = availFileKeys.intersect(data.keySet())
        if (foundFileKeys.size() == 0)
          exit 1, "TSV header must contain one of the file input keys: ${availFileKeys}}"
        if (foundFileKeys.size() > 1)
          exit 1, "TSV header cotains conflicting file keys: ${foundFileKeys}}"
        String fileKey = foundFileKeys.toList()[0]

        def file1 = data[fileKey]
        if (!file(file1).exists()) exit 1, "File does not exist: ${file1}"
        data.remove(fileKey)  // to be replaced with "fastq1" or "bam"

        if ((fileKey == 'fastq' || fileKey == 'fastq1') && !Utils.isFq(file1))
          exit 1, "fastq input file ${file1} does not have proper extension"
        if ((fileKey == 'bam') && !Utils.hasExtension(file1, '.bam'))
          exit 1, "bam input file ${file1} does not have proper extension"

        if (!data.containsKey('status')) data.status = 'N'  // normal sample by default
        data.status = Utils.returnStatus(data.status)

        if (file(file1).isDirectory()) {
          println("Path ${file1} is a directory, searching it for FastQ files")
          List fastqsFromDir = Utils.extractFastqFromDir(file1)
          fastqsFromDir.each { d ->
            data.each { d[it.key] = it.value }
            samples << d
          }
        } else {
          if (Utils.isFq(file1)) {
            samples << extractFastq(file1, data)
          }
          else if (hasExtension(file1, '.bam')) {
            samples << extractBam(file1, data)
          } else {
            "No recognisable extention for input file: ${file1}"
          }
        }
        if (!data.containsKey('sample')) data.sample = file1 -~/\.\w+$/
        if (!data.containsKey('patient')) data.patient = data.sample

        println "Parsed data: ${data}"
      }
    }

    samples
  }

  static Map extractBam(bam, data) {
    println("Path ${bam} is a BAM file")

    def bai = null
    if (data.containsKey('file2')) bai = data.file2
    if (data.containsKey('bai')) bai = data.bai
    if (bai) {
      if (!checkBai(bai, bam)) exit 1, "Second file ${bai} does not have the .bai extension"
      if (!bai.startsWith(bam -~/\.bam/) || !bai.startsWith(bam -~/\.BAM/))
        exit 1, "BAI file ${bai} does not match the BAM file ${bam}"
    }
    else // finding by extension
      bai = data.find { checkBai(it.value, bam) }

    if (bai && !file(bai).exists()) exit 1, "BAI file ${bai} does not exist"

    if (!bai) {
      bai = bam + ".bai"
      if (!file(bai).exists()) bai = (bam -~ /\.\w+$/) + ".bai"
      if (!file(bai).exists()) bai = null
    }

    data.bam = bam
    data.bai = bai
    return data
  }

  static Map extractFastq(fastq1, data) {
    println "Path ${fastq1} is a FastQ file"

    def fastq2 = null
    if (data.containsKey('file2')) fastq2 = data.file2
    if (data.containsKey('fastq2')) fastq2 = data.fastq2
    if (fastq2)
      if (!Utils.isFq(fastq2)) exit 1, "Second file input ${fastq2} is not FastQ"
    else  // finding by extension
      fastq2 = data.find { Utils.isFq(it.value) && it.value != fastq1 }

    if (!fastq2.exists()) error "Not found R2 FastQ file '${fastq2}'"

    if (fastq2) {
      println("Paths are FastQ files: ${fastq1} and ${fastq2}")
      if (!file(fastq2).exists()) exit 1, "Missing file in the TSV: ${fastq2}"
    } else {
      if (!fastq1.getName().contains('_R1')) exit 1, "Can't find R2 match for ${fastq1} " +
          "as the file name doesn't contain _R1. Workaround is to specify both R1 and R2 FastQ files " +
          "in the TSV with fastq1 and fastq2 keys."
      fastq2 = file(fastq1.toString().replace('_R1', '_R2'))
      if (!fastq2.exists()) error "Not found R2 FastQ file '${fastq2}'"
      println("Path ${fastq1} is a R1 fastq file. Found adjacent R2 file: ${fastq2}.")
    }

    data.fastq1 = fastq1
    if (fastq2) data.fastq2 = fastq2
    return data
  }

  static def checkTumorNormal(samples) {
    Map tumoursByPatient = [:]
    Map normalsByPatient = [:]

    samples.each {
      if (it.status == 0) normalsByPatient[it.patient] = it
      else                tumoursByPatient[it.patient] = it
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
  }

  static List extractFastqFromDir(pattern) {
    // Finding FASTQs from a directory pattern such as "my_samples/*/".
    // All FASTQ files in subdirectories are collected and emitted;
    // they must have _R1 and _R2 in their names.

    List samples = []

    List dirs = [file(pattern, type: "dir")].flatten()
    dirs.each {
      for (fastq1 in file("${it}/**_R1*.{fq,fastq}.gz")) {
        assert fastq1.getName().contains('_R1')
        def fastq2 = file(fastq1.toString().replace('_R1', '_R2'))
        if (!fastq2.exists()) error "Path '${fastq2}' not found"
        // the last name of the sampleDir is assumed to be a unique sample id
        String sampleId = fastq1.getParent().getName().toString()
        def (String flowcell, int lane) = Utils.flowcellLaneFromFastq(fastq1)
        GString rgId = "${flowcell}.${sampleId}.${lane}"
        samples << [sample: sampleId, rgId: rgId, fastq1: fastq1, fastq2: fastq2]
      }
    }
    samples
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

  static def startMessage(log, workflow, config, params, inputPath) {
    // Minimal information message
    log.info "Command line: " + workflow.commandLine
    log.info "Profile     : " + workflow.profile
    log.info "Project dir : " + workflow.projectDir
    log.info "Launch dir  : " + workflow.launchDir
    log.info "Work dir    : " + workflow.workDir
    log.info "Out dir     : " + params.outDir
    log.info "Input path  : " + inputPath
    log.info "Cont engine : " + workflow.containerEngine
    log.info "Executor    : " + config.process.executor
    log.info "Genome      : " + params.genome
    log.info "Genomes dir : " + params.genomes_base
  }

  static def endMessage(log, workflow, config, params, inputPath) {
    startMessage(log, workflow, config, params, inputPath)
    log.info "Completed at: " + workflow.complete
    log.info "Duration    : " + workflow.duration
    log.info "Success     : " + workflow.success
    log.info "Exit status : " + workflow.exitStatus
    log.info "Error report: " + (workflow.errorReport ?: '-')
  }
}