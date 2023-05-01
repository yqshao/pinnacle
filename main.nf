#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {entry as benchmark} from './nextflow/benchmark.nf'
include {entry as we_demo} from './nextflow/we_demo.nf'

workflow {we_demo()}
