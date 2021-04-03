# Teoroo Interactomic Potential Suite (TIPS)

TIPS is a set of tools for performing simulation with interatomic potentials.
TIPS is designed to streamline the composing of active learning workflows for
interatomic potentials.

## Installation

TIPS is yet not published 
A typical usage of TIPS is to use it for composing Nextflow workflows. In this
case, it's sufficient to include the tips module in your nextflow script, e.g.:

``` nextflow
nextflow.enable.dsl=2
include { pinnDs, pinnModel, aseSampler } from tips.nf
```

TIPS also includes a python library to to install it:

``` sh
pip install https://github.com/yqshao/tips.git
```

The python library also provides a command line utility for tasks like dataset
screening, conversion, etc. (see [here](https://tips.yqsho.github.io/latest/CLI)
for the documentation)

## Quick start

In a nutshell, TIPS modulize the ML problem as training processes (where models
are trained form datasets) and sampling processes (where datasets are generated
by models). In the simplist case, the training of a MLP can be represented as
such a workflow:

``` nextflow
nextflow.enable.dsl=2
include { aseDb, pinnTrainer, aseSampler } from tips.nf

workflow {
    data = aseDb('demo.db') 
    model = pinnTrainer(data, 'pinet.yml') 
    traj = aseSampler(model, 'npt.yml')
}
```

Here, the `activeLearn` function constructs an 

``` nextflow
nextflow.enable.dsl=2
include { dsConvert, pinnTrainer, aseSampler,  activeLearn } from tips.nf

workflow {
    data = dsConvert('demo.db', preprocess:true) 
    models, traj = activeLearn(
        pinnModel('pinet.yml'),
        cp2kModel('spce.lmp'), 
        aseSampler('npt.yml'),
        algo:'qbc', committee:3, fracTol:0.02, 
        errorTol:{energy: 0.01, force: 0.01})
}
```

The workflows are composed as an extention of the Nextflow workflow language, 
making the workflows extremely portable and easy to scale for a variety of 
computational resources. In the above case, it's trivial to run the workflow
on a HPC cluster with the below config file

``` nextflow
process {
    withlabel: 'pinn' {cpus = 20}
    withlabel: 'cp2k' {cpus = 40}
    executor='slurm'
}
singularity {
    enabled = true
}
```

While designed with active learning in mind, TIPS can also be used for writing
workflows involving *ab initio*, classical FFs or ML potentials.
