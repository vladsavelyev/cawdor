#!/usr/bin/env nextflow

if (params.help) exit 0, helpMessage()

genomeFasta = Utils.findRefFile(params, 'genomeFasta')
genomeIndex = Utils.findRefFile(params, 'genomeIndex')
genomeDict  = Utils.findRefFile(params, 'genomeDict')
purpleHet   = Utils.findRefFile(params, 'purpleHet')
purpleGC    = Utils.findRefFile(params, 'purpleGC')
dbsnp       = Utils.findRefFile(params, 'dbsnp')
intervals   = Utils.findRefFile(params, 'intervals')