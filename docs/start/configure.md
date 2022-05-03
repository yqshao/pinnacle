# Configure the workflow

## The nextflow config file

In your workflow folder, you will find a file named `nextflow.config` that
looks like this:

=== "Singularity"

    ```groovy
    profiles {
      standard {
        params {
          lmp_cmd = 'mpirun -np ${task.cpus} lmp_mpi'
          cp2k_cmd = 'mpirun -np ${task.cpus} cp2k.popt'}
        process {
          errorStrategy='ignore'
          withLabel: tips {container='yqshao/tips:tips-latest'}
          withLabel: pinn {container='yqshao/tips:pinn-latest'}
          withLabel: cp2k {container='yqshao/tips:cp2k-latest'}
          withLabel: utils {container='yqshao/tips:utils-latest'}
          withLabel: lammps {container='yqshao/tips:lammps-latest'}}
        executor {
          name = 'local'
          cpus = 4}}}

    singularity {
      enabled = true
      autoMounts = true}
    ```

=== "Slurm"

    ```groovy
    profiles {
      standard {
        params {
          lmp_cmd = 'mpirun -np ${task.cpus} lmp_mpi'
          cp2k_cmd = 'mpirun -np ${task.cpus} cp2k.popt'}
        process {
          errorStrategy='ignore'
          withLabel: tips {container='yqshao/tips:tips-latest'}
          withLabel: pinn {container='yqshao/tips:pinn-latest'}
          withLabel: cp2k {container='yqshao/tips:cp2k-latest'}
          withLabel: utils {container='yqshao/tips:utils-latest'}
          withLabel: lammps {container='yqshao/tips:lammps-latest'}}
        executor {
          name = 'local'
          cpus = 4}}}

    singularity {
      enabled = true
      autoMounts = true}
    ```

The file specifies details about your computational resources that is
independent of the workflow, e.g., the executable for your package, queuing
system, the resource you wish to use, etc.

## Adding multiple profiles

In the above example, the config is written to the "standard" profile. If you
would like to share the same workflow across different computational resources,
you can add additional profiles to the config, and switch them using the
`-profile` argument when running, i.e.:

=== "Command"

    ```bash
    nextflow run main.nf -profile my_profile
    ```

Check the nextflow documentation for some examples.

## General recommandation

When the project is generated with `tips wizard`, a minimal configuration file
is automatically generated as in the above example. The default provided by TIPS
might not fit your need. For instance, you might prefer the binaries compiled by
your local HPC cluster than the singularity images, or you might want to allocate
different resources to different types of calculations.

TIPS curates a list of [profiles]() for some of the computational resources we
have access to. Those can be good starting points for you to adapt for your own
need. You are also welcome to [contribute]() your config if think it will be
useful for others.
