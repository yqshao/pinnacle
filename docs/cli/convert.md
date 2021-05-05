# Convert

Convert datasets between different formats.

## Usage

```
$ tips convert dsfile1 [dsfile2 ...] [options]
```

## Options

| Option [shorthand] | Default     | Description                      |
|--------------------|-------------|----------------------------------|
| `--output [-o]`    | `'dataset'` | name of the output dataset       |
| `--format [-f]`    | `'auto'`    | format of input dataset          |
| `--oformat [-of]`  | `'pinn'`    | format of output dataset         |
| `--(no)shuffle`    | `False`     | shuffle dataset before splitting |
| `--seed`           | `0`         | random seed if shufffle is used  |
| `--units`          | `'real'`    | see [LAMMPS](lammps.md)          |
| `--emap`           | `''`        | see [LAMMPS](lammps.md)          |

The input `dsfile` will be concatenated when converting. The outputs will be
splitted if the output file is speficied as (e.g `train:8,test:2`). Note that
when `--splits` is used, different `dsfiles` will be splitted separately before
concatenting the result to the output files.
