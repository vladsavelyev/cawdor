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

Cawdor consists of several subworkflows. Currently supported map.nf to align reads and get alignment QC, and variants.nf to call somatic SNVs, indels, SVs, and CNVs.

Inputs can be either FastQ files or aligned BAM files.

To run, specify the input directory containing files (can be in subdirectories) with `--sampleDir`, or a TSV file with rows corresponding to samples with `--samples`. Also specify the output dir with `--outDir` and the genome build version with `--genome`

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


