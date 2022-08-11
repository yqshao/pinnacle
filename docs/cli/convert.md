# tips convert

Convert one or several datasets into another format. See
[tips.io](../python/io/#available-formats) for the formats available.

## Usage

```
pinn convert [options] ds1 ds2 ...
```

## Options

| Option [shorthand] | Default     | Description                                      |
|--------------------|-------------|--------------------------------------------------|
| `--output [-o]`    | `'dataset'` | name of the output dataset                       |
| `--fmt [-f]`       | `'auto'`    | format of input dataset                          |
| `--ofmt [-of]`     | `'runner'`  | format of output dataset                         |
| `--emap [-em]`     | `None`      | map the elements according to a LAMMPS data file |
