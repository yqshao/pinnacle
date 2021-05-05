# The Interactomic Potential Suite (TIPS)

TIPS is a set of tools for performing simulation with interatomic potentials.
TIPS is designed to streamline the composing of active learning workflows for
interatomic potentials.

**!!WARNING!!** this is a project in active development, nothing is guaranteed to work. 

## Quick start

Tips are contains a set of nextflow scripts that can be used without installation.

``` bash
nextflow run yqshao/tips --initDs 'test.xyz'
```

You can find how to install and extend TIPS in the
[documentation](https://yqshao.github.io/tips/workflow/overview/). TIPS also
provides a command line utility for tasks like dataset screening, conversion.
More information can be found here [here](https://yqshao.github.io/tips/cli/).

The workflows are composed in the Nextflow language, making the workflows
portable and easy to scale for a variety of computational resources. In the
above case, it's trivial to run the workflow on a HPC cluster with the below
config file

``` nextflow
process {
    withlabel: 'pinn' {cpus = 20}
    withlabel: 'lammps' {cpus = 1}
    executor='slurm'
}

singularity {
    enabled = true
}
```
