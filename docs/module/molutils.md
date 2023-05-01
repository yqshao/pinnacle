# Molecule Building Utilities

This module implements workflows related to the building of geometries or
building geometries as starting points of simulations.

## mol2box

The **tag** parameter should be a semicolon-separated list of SMILES strings,
with the number of molecules indicated by a number, separated with comma. For
example, `o,64` means a box with 64 water molecule, and `O,32;COO,32` means an
equimolar mixture of water and ethonol.

The **box** parameter can be a single float for a cubic box, a
semicolon-separated list for an orthogonal cell.

The **seed** parameter specified the random number generator used in packmol.

The 3D conformation will be build and hydrogen atoms added by [openbabel]. The
molecules will then be packed into the specified box with [packmol].

[openbabel]: https://openbabel.org/docs/dev/Command-line_tools/babel.html#options
[packmol]: https://m3g.github.io/packmol/download.shtml

### Input channels

| Channel | Type   | i/o[idx]     | Note                             |
|---------|--------|--------------|----------------------------------|
| `name`  | `val`  | `in[0]`      | an id to identify the process    |
| `tag`   | `val`  | `in[1]`      | a string specifying the geometry |
| `box`   | `val`  | `in[2]`      | a string specifying the box size |
| `seed`  | `val`  | `in[3]`      | random number generator seed     |
| `name`  | `val`  | `out.*[0]`   | same as input                    |
| `geo`   | `file` | `out.geo[1]` | the geometry packed with packmol |
| `log`   | `file` | `out.log[1]` | packmol input and logs           |

??? "Source code"

    ```groovy
    --8<-- "nextflow/module/molutils.nf"
    ```

