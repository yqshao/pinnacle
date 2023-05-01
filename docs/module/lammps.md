# LAMMPS workflows

## lammpsMD

The `lammpsMD` process runs a LAMMPS simulation is specified with an input file,
and the set of auxiliary files. The process takes the `params.lmp_cmd` to
specify the lammps binary, which can be specified during runtime or in
`nextflow.config`.

### Channels

| Channel | Type | i/o[idx] | Note                                      |
|---------|------|----------|-------------------------------------------|
| name    | val  | `in[0]`  | an id to identify the process             |
| input   | file | `in[1]`  | LAMMPS input file                         |
| aux     | file | `in[2]`  | auxiliary files (force field, data, etc.) |
| name    | val  | `out[0]` | same as output                            |
| traj    | file | `out[1]` | trajectory in `.dump` format              |
| log     | file | `out[2]` | LAMMPS log                                |
| restart | file | `out[3]` | LAMMPS restart files in `.restart`        |

??? "Source code"

    ```groovy
    --8<-- "nextflow/module/lammps.nf"
    ```
