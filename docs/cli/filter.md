# filter 

Filter the dataset according to certain criterion.

## Usage

```
$ tips filter dsfile1 [dsfile2 ...] [options]
```

## Description

| Name [shorthand]    | Default     | Description                        |
|---------------------|-------------|------------------------------------|
| `--output [-o]`     | `'dataset'` | name of the output dataset         |
| `--format [-f]`     | `'auto'`    | format of input dataset            |
| `--oformat [-of]`   | `'pinn'`    | format of output dataset           |
| `--algo [-a]`       | `'naive'`   | Algorithm to perform the filtering |
| `--val-max [-vmax]` | `''`        | Maximum value, see below           |
| `--val-min [-vmin]` | `''`        | Minimum value, see below           |
| `--abs-max [-amax]` | `''`        | Maximum absolute value, see below  |
| `--abs-min [-amin]` | `''`        | Minimum absolute value, see below  |
| `--units`           | `'real'`    | see [LAMMPS](lammps.md)            |
| `--emap`            | `''`        | see [LAMMPS](lammps.md)            |
    
The filter tags are interpreted as comma separated `key:val` pairs. Structure is
excluded if **any** of the condition is met, for instance `--vmax "e:-20",
--amax "f:10"` means any structure with energy higher than `-20` or force
component larger than `10` will be excluded. 

In some cases is might be useful to alter this behavior, and get instead the
structures with at least one force component larger than 10. This can be done 
`--amax "!f:10"`.

## Algorithms

- **naive**: in the simplest case the filter is applied to the values in the dataset.
- **qbc**: the qbc filter takes the input files as copies of the same structure
  with different labels. It will check that all the datasets have the same
  structure and filter according to the standard deviation of the labels. (when
  only two datasets are used, this is equalalent as filtering according to the
  RMSE).
