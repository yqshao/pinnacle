#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.datasets = './datasets/*.data'
params.pinn_inp = './inputs/pinet.yml'
params.ase_init = './inputs/init.xyz'

// below are addition parameters that might be ajusted
params.repeats = 5
params.pinn_flags = '--train-steps 5000000 --log-every 5000 --ckpt-every 100000 --batch 1 --max-ckpts 1 --shuffle-buffer 1000 --init'
params.ase_flags = '--ensemble npt --T 373 --t 100 --dt 0.5 --log-every 10'
params.rdf_flags = '--tags O-O,O-H --rc 5 '
params.log_flags = '--tags density'

include { pinnTrain } from './pinn.nf' addParams(publish: 'models')
include { aseMD } from './ase.nf' addParams(publish: 'trajs')
include { rdf; mdlog } from './analysis.nf' addParams(publish: 'analyses')

workflow {
  channel.fromPath(params.datasets) \
    | combine(channel.fromPath(params.pinn_inp)) \
    | combine(channel.of(1..params.repeats)) \
    | map { ds, inp, seed -> ["$ds.baseName-$inp.baseName-$seed", ds, inp, "--seed $seed $params.pinn_flags"] } \
    | pinnTrain \
    | map { name, model -> [name, model, file(params.ase_init), params.ase_flags] } \
    | aseMD

  aseMD.out | map {name, traj -> [name, traj, params.rdf_flags]} | rdf
  aseMD.out | map {name, traj -> [name, traj, params.log_flags]} | mdlog
}
