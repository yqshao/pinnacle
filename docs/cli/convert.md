# tips convert

Convert one or several datasets into another format. See
[tips.io](../python/io/#available-formats) for the formats available.

## Usage

```
tips convert [options] ds1 ds2 ...
```

## Options

| Option [shorthand] | Default    | Description                                      |
|--------------------|------------|--------------------------------------------------|
| `--fmt [-f]`       | `'auto'`   | format of input dataset                          |
| `--emap [-em]`     | `None`     | map the elements according to a LAMMPS data file |
| `--output [-o]`    | `'output'` | name of the output dataset                       |
| `--ofmt [-of]`     | `'extxyz`  | format of output dataset                         |
| `--shuffle`        | `False`    | shuffle the dataset (only for indexable formats) |
| `--seed`           | `0`        | random seed for shuffling                        |
