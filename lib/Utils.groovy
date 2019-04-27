import static nextflow.Nextflow.file
import nextflow.Channel

class Utils {
  
  static def defineReferenceMap(params) {
    if (!(params.genome in params.genomes)) exit 1, "Genome ${params.genome} not found in configuration"
    return [
      // genome reference dictionary
      'genomeDict'       : Utils.checkParamReturnFile(params, "genomeDict"),
      // FASTA genome reference
      'genomeFasta'      : Utils.checkParamReturnFile(params, "genomeFasta"),
      // genome .fai file
      'genomeIndex'      : Utils.checkParamReturnFile(params, "genomeIndex"),
      // BWA index files
      'bwaIndex'         : Utils.checkParamReturnFile(params, "bwaIndex"),
      // intervals file for spread-and-gather processes
      'intervals'        : Utils.checkParamReturnFile(params, "intervals"),
      // for variant annotation
      'dbsnp'            : Utils.checkParamReturnFile(params, "dbsnp"),
      'dbsnpIndex'       : Utils.checkParamReturnFile(params, "dbsnpIndex")
    ]
  }

  static def checkExactlyOne(list) {
    def n = 0
    list.each{n += it ? 1 : 0}
    return n == 1
  }

  static def checkParamReturnFile(params, item) {
    params."${item}" = params.genomes[params.genome]."${item}"
    if (!params."${item}") exit 1, "params do not have item ${item}, params=${params}"
    return file(params."${item}")
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
  static def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
  }

  // Check file extension
  static def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
  }

  // Check if the file extention corresponds to fastq
  static def isFq(it) {
    hasExtension(it, 'fastq.gz') || hasExtension(it, 'fq.gz')
  }

  // Check parameter existence
  static def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
      println("Unknown parameter: ${it}")
      return false
    }
    return true
  }

  // Compare each parameter with a list of parameters
  static def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
  }  

  // Loop through all the references files to check their existence
  static def checkReferenceMap(referenceMap) {
    referenceMap.every {
      referenceFile, fileToCheck ->
      Utils.checkRefExistence(referenceFile, fileToCheck)
    }
  }  

  // Loop through all the references files to check their existence
  static def checkRefExistence(referenceFile, fileToCheck) {
    if (fileToCheck instanceof List) return fileToCheck.every{ Utils.checkRefExistence(referenceFile, it) }
    def f = file(fileToCheck)
    // this is an expanded wildcard: we can assume all files exist
    if (f instanceof List && f.size() > 0) return true
    else if (!f.exists()) {
      println  "Missing references: ${referenceFile} ${fileToCheck}"
      return false
    }
    return true
  }

  // Channeling the TSV file containing BAM.
  // Format is: "subject status sample bam bai"
  static def extractBams(tsvFile) {
    Channel.from(tsvFile)
      .splitCsv(sep: '\t')
      .map { row ->
        SarekUtils.checkNumberOfItem(row, 6)
        def idPatient = row[0]
        def status    = Utils.returnStatus(row[1].toInteger())
        def idSample  = row[2]
        def bamFile   = Utils.returnFile(row[3])
        def baiFile   = Utils.returnFile(row[4])

        if (!SarekUtils.hasExtension(bamFile, ".bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
        if (!SarekUtils.hasExtension(baiFile, ".bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

        return [ idPatient, status, idSample, bamFile, baiFile ]
      }
  }
}

