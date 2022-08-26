# tips utils

This commnad collects several tools useful in preparing simulation input files.

## mkcp2kinp

This command takes a CP2K `inp` file and inserts the structural data to it. It
supprts single-shot mode, where one output `cp2k.inp` will be made; and a
subsample mode, where subsampled dataset will be used to generate several input
files.

### Usage
```
tips utils mkcp2kinp [options] inp dataset
```

### Options

| Option [shorthand] | Default       | Description                                       |
|--------------------|---------------|---------------------------------------------------|
| `--fmt [-f]`       | `'auto'`      | format of input dataset                           |
| `--emap [-em]`     | `None`        | map the elements according to a LAMMPS data file  |
| `--idx`            | `-1`          | index of the used datum, used in single-shot mode |
| `--subsample`      | `False`       | activates the subsample mode                      |
| `--strategy`       | `'uniform'`   | one of 'uniform' or 'sorted'                      |
| `--nsample`        | `None`        | number to subsample                               |
| `--psample`        | `None`        | percentage to subsample                           |
| `--sort-Key`       | `'force_std'` | key used in the sorted scheme                     |

