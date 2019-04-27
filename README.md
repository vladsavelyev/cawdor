# Cawdor

**Cancer analysis workflow for DNAseq or RNAseq**

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

### Inputs

Input files can be either raw FastQ files or aligned BAM files. To specify input files, you can either provide a directory, or a TSV file.

#### Samples directory

To run, you can specify the input directory with `--sampleDir`. The directory searched recursively for FastQ files that are named `*_R1_*.fastq.gz`, and a matching pair with `_R2_` instead of `_R1_`):

```bash
nextflow run align.nf --samplesDir /samples
```

For multiple patiens, organize the foloder into one subfolder for every sample:

```
ID
+--sample1
+------sample1_lib_flowcell-index_lane_R1_1000.fastq.gz
+------sample1_lib_flowcell-index_lane_R2_1000.fastq.gz
+------sample1_lib_flowcell-index_lane_R1_1000.fastq.gz
+------sample1_lib_flowcell-index_lane_R2_1000.fastq.gz
+--sample2
+------sample2_lib_flowcell-index_lane_R1_1000.fastq.gz
+------sample2_lib_flowcell-index_lane_R2_1000.fastq.gz
+--sample3
+------sample3_lib_flowcell-index_lane_R1_1000.fastq.gz
+------sample3_lib_flowcell-index_lane_R2_1000.fastq.gz
+------sample3_lib_flowcell-index_lane_R1_1000.fastq.gz
+------sample3_lib_flowcell-index_lane_R2_1000.fastq.gz
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

#### Samples TSV file

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

#### Other options

Also set the output dir with `--outDir`, and the genome build version with `--genome`.

```
nextflow run map.nf --sampleDir ../Sarek/Sarek-data/testdata/tiny --outDir Results --genome smallGRCh37
```


To run on a specific infrastructure like NCI Raijin or Spartan, use a corresponding profile: 

```
nextflow run map.nf --sampleDir ../Sarek/Sarek-data/testdata/tin --outDir Results --genome smallGRCh37 -profile raijin
```

To run on NCI on a local node using 1 cpu, you can use `-process.*` nextflow options:

```
nextflow run map.nf --sampleDir ../Sarek/Sarek-data/testdata/tin --outDir Results --genome smallGRCh37 -profile raijin -process.cpus=1 -process.executor=local
```


