# Cawdor

**Cancer analysis workflow for (WGS) DNAseq or RNAseq**

[![Build Status](https://travis-ci.com/vladsaveliev/cawdor.svg?branch=master)](https://travis-ci.com/vladsaveliev/cawdor)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)


## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a portable manner. The structure is inspired by [bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen), a python-based NGS analysis framework; the implementation is inspired by [Sarek](https://github.com/SciLifeLab/Sarek), a nextflow-based cancer analysis workflow; some approaches are borrowed from [Hartwig Medical Foundation pipeline](https://github.com/hartwigmedical/hmftools/).


## Installation

For installation, follow the nextflow and nf-core documentation:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Troubleshooting](https://nf-co.re/usage/troubleshooting)


## Usage

Cawdor consists of several subworkflows: `align.nf` to align reads and get alignment QC, `somatic.nf` to call somatic variants (SNVs, indels, SVs, and CNVs), `germline.nf` to call germline variants, and `postprocess.nf` to annotate and prioritise variants, generate reports and QC.

The typical command for running the pipeline is as follows:

```bash
nextflow run align.sh --samplesDir /samples -profile raijin --outDir Results --genome GRCh37
```

This will launch the pipeline with the `raijin` cluster configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work/            # Directory containing the nextflow working files
Results/         # Finished results (configurable, see below)
.nextflow.log    # Log file from Nextflow
.nextflow/       # Folder with other Nextflow hidden files
```


## Specifying input data

Input files can be either raw FastQ files or aligned BAM files. To specify input files, you can either provide a directory, or a TSV file.

### Samples directory

To run, you can specify the input directory with `--sampleDir`. The directory searched recursively for FastQ files that are named `*_R1_*.fastq.gz`, and a matching pair with `_R2_` instead of `_R1_`):

```bash
nextflow run align.nf --samplesDir /samples
```

For multiple patients, organize the folder into one subfolder for every sample:

```
ID
+--sample1
+-----sample1_lib_flowcell-index_lane_R1_1000.fastq.gz
+-----sample1_lib_flowcell-index_lane_R2_1000.fastq.gz
+-----sample1_lib_flowcell-index_lane_R1_1000.fastq.gz
+-----sample1_lib_flowcell-index_lane_R2_1000.fastq.gz
+--sample2
+-----sample2_lib_flowcell-index_lane_R1_1000.fastq.gz
+-----sample2_lib_flowcell-index_lane_R2_1000.fastq.gz
+--sample3
+-----sample3_lib_flowcell-index_lane_R1_1000.fastq.gz
+-----sample3_lib_flowcell-index_lane_R2_1000.fastq.gz
+-----sample3_lib_flowcell-index_lane_R1_1000.fastq.gz
+-----sample3_lib_flowcell-index_lane_R2_1000.fastq.gz
```

FastQ filename structure:

- `sample_lib_flowcell-index_lane_R1_1000.fastq.gz` and
- `sample_lib_flowcell-index_lane_R2_1000.fastq.gz`

Where:

- `sample` = sample id
- `lib` = indentifier of libaray preparation
- `flowcell` = identifyer of flow cell for the sequencing run
- `lane` = identifier of the lane of the sequencing run

Read group information will be parsed from fastq file names according to this:

- `RGID` = "sample_lib_flowcell_index_lane"
- `RGPL` = "Illumina"
- `PU` = sample
- `RGLB` = lib

### Samples TSV file

Another option is to specify a TSV file with rows corresponding to samples with `--samples`.

```bash
nextflow run align.nf --samples samples.tsv
```

The TSV file should have at least one tab-separated line:

```
SUBJECT_ID_1	0	SAMPLE_1_N	1	/samples/normal1_1.fastq.gz	/samples/normal1_2.fastq.gz
SUBJECT_ID_1	1	SAMPLE_1_T	3	/samples/tumor1_1.fastq.gz	/samples/tumor1_2.fastq.gz
SUBJECT_ID_2	0	SAMPLE_2_N	2	/samples/normal2_1.fastq.gz	/samples/normal2_2.fastq.gz
SUBJECT_ID_2	1	SAMPLE_2_T	4	/samples/tumor2_1.fastq.gz	/samples/tumor2_2.fastq.gz
```

The columns are:

1. Subject (batch) id
2. Status: 0 if normal, 1 if tumor
3. Sample id: actual text representation of the type of the sample
4. Lane ID - used when the sample is multiplexed on several lanes
5. First set of reads
6. Second set of reads

To run from BAM file, create a 5-column TSV file:

```
SUBJECT_ID_1	0	SAMPLE_1_N	1	/samples/normal_1.bam
SUBJECT_ID_1	1	SAMPLE_1_T	3	/samples/tumor_1.bam
SUBJECT_ID_2	0	SAMPLE_2_N	2	/samples/normal_2.bam
SUBJECT_ID_2	1	SAMPLE_2_T	4	/samples/tumor_2.bam
```

Another option is to specify one folder for each sample. Useful when you have many lanes:

```
SUBJECT_ID_1	0	SAMPLE_1_N	1	/samples/sample_n_1
SUBJECT_ID_1	1	SAMPLE_1_T	3	/samples/sample_t_1
SUBJECT_ID_2	0	SAMPLE_2_N	2	/samples/sample_n_2
SUBJECT_ID_2	1	SAMPLE_2_T	4	/samples/sample_t_2
```

## Reference genomes

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

### `--genome` (using iGenomes)

There are 31 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config). Common genomes that are supported are:

* Human
  * `--genome GRCh37`
  * `--genome hg38`

Presets exist for `raijin` and `spartan` environments. For other machines, provide the location of the genomes with `--genomes_dir` option, with the directory having the following structure:

```
bwaIndex         = "${params.genome_base}/${params.genome}/${params.genome}.fa.{amb,ann,bwt,pac,sa}"
genomeDict       = "${params.genome_base}/${params.genome}/${params.genome}.dict"
genomeFasta      = "${params.genome_base}/${params.genome}/${params.genome}.fa"
genomeIndex      = "${params.genome_base}/${params.genome}/${params.genome}.fa.fai"
intervals        = "${params.genome_base}/${params.genome}/wgs_calling_regions_CAW.list"
dbsnp            = "${params.genome_base}/${params.genome}/dbsnp-151.vcf.gz"
dbsnpIndex       = "${params.genome_base}/${params.genome}/dbsnp-151.vcf.gz.tbi"
vepCacheVersion  = "94"
```

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `raijin`
  * Uses NCI PBSPro scheduler as exeutor, also knows about available resourses, default location to the conda environment and reference genomes. 
* `spartan`
  * Uses Spartan Slurm scheduler as exeutor, also knows about available resourses, default location to the conda environment and reference genomes.   
* `awsbatch`
  * A generic configuration profile to be used with AWS Batch.
* `conda`
  * A generic configuration profile to be used with [conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`vladsaveliev/cawdor`](http://hub.docker.com/r/vladsaveliev/cawdor/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`vladsaveliev/cawdor`](http://hub.docker.com/r/vladsaveliev/cawdor/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

Example: running on NCI Raijin: 

```
nextflow run align.nf --sampleDir ../Sarek/Sarek-data/testdata/tin --outDir Results --genome smallGRCh37 -profile raijin
```

You can also run on NCI on a local node, using 1 cpu - just set up `-process.*` nextflow options:

```
nextflow run align.nf --sampleDir ../Sarek/Sarek-data/testdata/tin --outDir Results --genome smallGRCh37 -profile raijin -process.cpus=1 -process.executor=local
```


## Other command line parameters

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.




