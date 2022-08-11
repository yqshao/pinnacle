# CP2K workflows

## cp2k

The CP2K process runs, a CP2K MD is specified with an input file, and the set of
auxiliary files. For conviniece there's also workflows like [cp2kMD](#cp2kMD)
that takes extra channels and prepares the input files. Those processes shares
two parameters that can be specified during runtime or in `nextflow.config`:

- `params.cp2k_cmd`: the command for invoking cp2k
- `params.cp2k_aux`: the path (with wildcards) to the potential, basis, ...

### Channels

| Channel       | Type | Note                                      |
| ------------- | ---- | ----------------------------------------- |
| (in) name     | val  | an id to identify the process             |
| (in) input    | file | CP2k input file                           |
| (in) aux      | file | Auxiliary files                           |
| (out) traj    | file | Trajectory (pos, frc, ener, cell, stress) |
| (out) log     | file | CP2K log                                  |
| (out) restart | file | CP2K restart file                         |

## cp2kMD

This is a shortcut to run an CP2K MD form given initial geometry, by taking an
input file. The file will be inserted to the input as the initial gemoetry. The
geometry must be recognizable by `tips.io`. In case that a multi-frame trajectry
is used, the last frame will be used.

### Channels

| Channel[^out] | Type | Note                                         |
| ------------- | ---- | -------------------------------------------- |
| (in) name     | val  | an id to identify the process                |
| (in) input    | val  | a CP2k input file                            |
| (in) init     | file | a geometry file, to be inserted to the input |

[^out]: Output channels are the same as [cp2k](#cp2k).

??? "Source code"

    ```groovy
    --8<-- "nextflow/cp2k.nf"
    ```
