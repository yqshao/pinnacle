#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// This is an exmple of benchmarking script for a liquid system, multiple
// datasets and models may be specied but the same analysis is performed
// for each trained model.

// to prepare a dataset from a RuNNer-formatted dataset:
// ``
// to prepare a dataset from a DeePMD-formatted dataset:
//

// the paramsters include those files:
// dataset/*.data, several input datasets
// inputs/*.yml, several model parameters
// inputs/init.xyz, initial geometry for
params.datasets = './datasets/*.data'
params.pinn_inp = './inputs/pinet.yml'
params.asemd_init = './inputs/init.xyz'

// below are addition parameters that might be ajusted
params.repeats = 5
params.pinn_flags = '--train-steps 5000000 --log-every 5000 --ckpt-every 100000 --batch 1 --max-ckpts 1 --shuffle-buffer 1000 --init'
params.ase_flags = '--ensemble npt --T 373 --t 100 --dt 0.5 --every 10'
params.rdf_flags = '--tags O-O,O-H --rc 5 '
params.log_flags = '--tags density'

include { pinnTrain } from './pinn.nf'
include { aseMD } from './ase.nf'
include { rdf; mdlog } from './analysis.nf'

workflow {
  // combine
  channel.fromPath(params.datasets)
    .combine(channel.fromPath(params.pinn_inp))
    .combine(channel.of(1..params.repeats))
    .multiMap{ ds, inp, seed ->
      ds: ds
      inp: inp
      flag: "--seed $seed $params.pinn_flags"
      name: "$ds.baseName-$inp.baseName-$seed"
      path: "models/$ds.baseName-$inp.baseName-$seed"
    }
    .set {ch}

  pinnTrain(ch.name, ch.ds, ch.inp, ch.flag, ch.path).multiMap{
    name, model ->
    name: name
    model: model
    path: "trajs/$name"
  }.set {ch}

  init = file(params.asemd_init)
  flags = params.ase_flags
  aseMD(ch.name, ch.model, init, flags, ch.path).traj.multiMap{
    name, traj ->
    name: name
    traj: traj
    path: "analyses/$name"
  }.set {ch}

  rdf(ch.name, ch.traj, params.rdf_flags, ch.path)
  mdlog(ch.name, ch.traj, params.log_flags, ch.path)
}
