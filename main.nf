#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {wf_benchmark} from './nextflow/benchmark.nf'

workflow benchmark {wf_benchmark()}
