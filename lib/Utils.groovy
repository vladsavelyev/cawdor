import static nextflow.Nextflow.file
import nextflow.Channel
import static nextflow.Nextflow.exit

class Utils {
  
  static def findRefFile(params, key) {
    def genomes = params.genomes
    def genome = params.genome

    if (!(genome in genomes))
      exit 1, "Genome ${genome} not found in configuration. Available: ${genomes}"

    def path = genomes[genome]."${key}"
    if (!path) {
      path = genomes["default"]."${key}"
      if (!path) exit 1, "Missing key ${key} for genome ${genome}"
    }

    path = file(path)
    if (!Utils.checkRefExistence(key, path)) {
      println "Not found reference file ${key}: ${path}"
      return null
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
      println "Missing reference files: ${key} ${fileToCheck}"
      return false
    }
    return true
  }

  static boolean checkExactlyOne(list) {
    def n = 0
    list.each{n += it ? 1 : 0}
    return n == 1
  }

  // Return status [0,1]
  // 0 == Normal, 1 == Tumor
  static def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
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
    Channel.from(tsvFile)
      .splitCsv(sep: '\t')
      .map { row ->
        Utils.checkNumberOfItem(row, 4)
        def idPatient = row[0]
        def status    = Utils.returnStatus(row[1].toInteger())
        def idSample  = row[2]
        def bamFile   = Utils.returnFile(row[3])
        def baiFile   = bamFile + ".bai"

        if (!Utils.hasExtension(bamFile, ".bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
        if (!Utils.hasExtension(baiFile, ".bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

        return [ idPatient, status, idSample, bamFile, baiFile ]
      }
  }
}

