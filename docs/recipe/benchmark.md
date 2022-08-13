# Benchmark workflow

Benchmarking is a typical usage of the TIPS framework. It works with a fixed
dataset and test the result of different simulations. The possible processes are
listed below, along with an annoted sample workflow.

=== "Flowchart"

    ```mermaid
    graph LR
    ds[dataset] --> train([Train])
    inp[input] --> train
    seeds --> train
    train --> model[Model]
    model --> md([MD])
    model --> vscan([Volum Scan])
    md --> trajs[Trajectory]
    trajs --> rdf([RDF])
    trajs --> mdlog([MDLog])
    ```

=== "main.nf"

    ```groovy
    // below are addition parameters that might be ajusted
    params.datasets = './datasets/*.data' // (1)
    params.pinn_inp = './inputs/pinet.yml'
    params.ase_init = './inputs/init.xyz'
    params.repeats = 5
    params.pinn_flags = '--train-steps 5000000'
    params.ase_flags = '--ensemble npt --T 373 --t 100 --dt 0.5 --log-every 10'
    params.rdf_flags = '--tags O-O,O-H --rc 5 '
    params.log_flags = '--tags density'

    include { pinnTrain } from './pinn.nf' addParams(publish: 'models')
    include { aseMD } from './ase.nf' addParams(publish: 'trajs')
    include { rdf; mdlog } from './analysis.nf' addParams(publish: 'analyses')

    workflow { // (2)
      channel.fromPath(params.datasets)                       \
        | combine(channel.fromPath(params.pinn_inp))          \
        | combine(channel.of(1..params.repeats))              \
        | map { ds, inp, seed ->                              \
                ["$ds.baseName-$inp.baseName-$seed", ds, inp, \
                 "--seed $seed $params.pinn_flags"] }         \
        | pinnTrain                                           \
        | combine(channel.fromPath(params.ase_init))          \
        | map { name, model, init ->                          \
                [name, model, init, params.ase_flags] }       \
        | aseMD // (3)

      aseMD.out | map {name, traj -> [name, traj, params.rdf_flags]} | rdf
      aseMD.out | map {name, traj -> [name, traj, params.log_flags]} | mdlog
    }
    ```

    1. The files can be specified with wildcards to match multiple files.
    2. The [`combine` operator](https://www.nextflow.io/docs/latest/operator.html#combine)
       builds `tuple`s with all possible combinations. Here, we want to benchmark all dataset
       with all model inputs, and run `repeats` for each model.
    3. to apply multiple analyses to the same trajectory, you can use the `.out` of the
       processes.

## Other links

- [List of available datasets]()
- [Sample model input files]()
- [Adding a custom analysis]()
