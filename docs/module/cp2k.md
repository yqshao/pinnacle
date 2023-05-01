# CP2K workflows

## Parameters

This module uses the following paramters:

| Parameter  | Default     | Note                    |
|------------|-------------|-------------------------|
| `publish`  | `cp2k`      | default publish folder  |
| `cp2k_cmd` | `cp2k.psmp` | command to run CP2K    |
| `cp2k_aux` | `null`      | path to auxiliary files |

## cp2k

The `cp2k` process is specified with an input file, and the set of auxiliary
files. 

### Channels

| Channel   | i/o[idx]         | Type | Note                                      |
|-----------|------------------|------|-------------------------------------------|
| `name`    | `in[0] `         | val  | an id to identify the process             |
| `input`   | `in[1] `         | file | CP2k input file                           |
| `aux`     | `in[2] `         | file | auxiliary files                           |
| `name`    | `out.*[0]`       | file | same as input                             |
| `traj`    | `out.traj[1]`    | file | trajectory (pos, frc, ener, cell, stress) |
| `log`     | `out.log[1]`     | file | CP2K log                                  |
| `restart` | `out.restart[1]` | file | CP2K restart file                         |

## cp2kMD

This is a proxy to run `cp2k` from a given initial geometry, by taking an input
file. The file will be inserted to the input as the initial gemoetry. The
geometry must be recognizable by `tips.io`. In case that a multi-frame trajectry
is used, the last frame will be used.

### Channels

| Channel   | i/o[idx]         | Type | Note                                      |
|-----------|------------------|------|-------------------------------------------|
| `name`    | `in[0] `         | val  | an id to identify the process             |
| `input`   | `in[1] `         | file | CP2k input file                           |
| `aux`     | `in[2] `         | file | auxiliary files                           |
| `name`    | `out.*[0]`       | file | same as input                             |
| `traj`    | `out.traj[1]`    | file | trajectory (pos, frc, ener, cell, stress) |
| `log`     | `out.log[1]`     | file | CP2K log                                  |
| `restart` | `out.restart[1]` | file | CP2K restart file                         |

??? "Source code"

    ```groovy
    --8<-- "nextflow/module/cp2k.nf"
    ```
