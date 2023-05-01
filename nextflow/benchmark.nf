#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.proj = 'benchmark'
params.ds = 'datasets/*'
params.inp = 'inputs/*'
params.flags = '--train-steps 3e6 --batch 30 --shuffle-buffer 3000'
params.seeds = 3

include {train} from './module/pinn.nf' addParams(publish: "$params.proj/models")

workflow entry {
  ch_ds = Channel.fromFilePairs("${params.ds}.{yml,tfr}")
  ch_inp = Channel.fromPath("${params.inp}.yml")
  ch_seed = Channel.of(1..params.seeds.toInteger())

  ch_ds \
    | combine(ch_inp) \
    | combine(ch_seed) \
    | map {ds, ds_files, inp, seed -> \
           ["$ds-${inp.baseName}-$seed", ds_files, inp, "$params.flags --seed $seed"]} \
    | train
}

// default entrypoint
workflow {entry()}
