# tips/lammps

This module implements the `sampler` and `labeller` process with
[LAMMPS](https://lammps.sandia.gov/) molecular dynamics code. LAMMPS allow for a
flexible input format. Typically, the binary is run with an input script
(`input.lmp` here) that holds all the information required for a simulation:

```shell
mpirun lmp_mpi -in input.lmp
```

while extra input files may be included in the script, e.g.:

```shell
# Include the force field, and geometry files
include input.init
read_data input.data
include input.setting
# The md setup (ensemble, etc) goes below
...
```

For both `lammpsSample` and `lammpsLabel`, the input script is required as
`inp`. Some assumption is made about the input script, following the convention
of moltemplate. The input script should read three files in order, specified by
`lmpInit`, `lmpData` and `lmpSetting`. They are saved as `input.init`,
`input.data` and `input.setting` respectively in the process.

The `init` file will be saved as `input.coord` for `lammpsSample` workflow and
 `ds` file as `input.dump` for `lammpsLabel`. It is assumed that the LAMMPS
 input script read those files and produces two files `output.log` and
 `output.dump`, which holds the thermo logs from lammps and the trajectory,
 respectively. In addition, the `emap` and `units` may be supplied to inform
 TIPS about the units and element specification in the LAMMPS simulations.

You can find some the example input files for LAMMPS in the [example
folder](https://github.com/yqshao/tips/tree/master/examples/explore-nacl).

## lammpsSample

| Inputs     | Default  | Description             |
|------------|----------|-------------------------|
| inp        | `null`   | LAMMPS input script     |
| init       | `null`   | initial structure       |
| lmpEmap    | `''`     | relabel LAMMPS elements |
| lmpUits    | `'real'` | LAMMPS units            |
| lmpInit    | `null`   | LAMMPS initialization   |
| lmpData    | `null`   | LAMMPS data             |
| lmpSetting | `null`   | LAMMPS setting          |

Outputs the sampled trajectory.

## lammpsLabel
Label a dataset with Lammps

| Inputs     | Default  | Description             |
|------------|----------|-------------------------|
| ds         | `null`   | dataset to label        |
| inp        | `null`   | LAMMPS input script     |
| lmpEmap    | `''`     | relabel LAMMPS elements |
| lmpUnits   | `'real'` | LAMMPS units            |
| lmpInit    | `null`   | LAMMPS initialization   |
| lmpData    | `null`   | LAMMPS data             |
| lmpSetting | `null`   | LAMMPS setting          |

Outputs the labelled dataset.
