# tips subsample

Subsample a dataset to get a subset

## Usage 

```
tips subsample [options] dataset
```

## options

| Option [shorthand] | Default       | Description                                      |
|--------------------|---------------|--------------------------------------------------|
| `--fmt [-f]`       | `'auto'`      | format of input dataset                          |
| `--emap [-em]`     | `None`        | map the elements according to a LAMMPS data file |
| `--output [-o]`    | `'output'`    | name of the output dataset                       |
| `--ofmt [-of]`     | `'extxyz`     | format of output dataset                         |
| `--strategy`       | `'uniform'`   | one of 'uniform' or 'sorted'                     |
| `--nsample`        | `None`        | number to subsample                              |
| `--psample`        | `None`        | percentage to subsample                          |
| `--sort-Key`       | `'force_std'` | key used in the sorted scheme                    |
