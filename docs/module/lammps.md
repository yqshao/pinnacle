# LAMMPS workflows

## lammps

The CP2K process runs, a CP2K MD is specified with an input file, and the set of
auxiliary files. The process takes the `params.lmp_cmd` to specify the lammps
binary, which can be specified during runtime or in `nextflow.config`.

### Channels

| Channel       | Type | Note                                      |
| ------------- | ---- | ----------------------------------------- |
| (in) name     | val  | an id to identify the process             |
| (in) input    | file | LAMMPS input file                         |
| (in) aux      | file | Auxiliary files (force field, data, etc.) |
| (in) publish  | val  | path to publish the output traj. and log  |
| (out) traj    | file | Trajectory in `.dump` format              |
| (out) log     | file | LAMMPS log                                |
| (out) restart | file | LAMMPS restart files in `.restart`        |

??? "Source code"

    ```groovy
    --8<-- "nextflow/lammps.nf"
    ```
