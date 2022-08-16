# tips convert

Convert one or several datasets into another format. See
[tips.io](../python/io/#available-formats) for the formats available.

## Usage

```
pinn convert [options] ds1 ds2 ...
```

## Options

| Option [shorthand] | Default       | Description                                      |
| ------------------ | ------------- | ------------------------------------------------ |
| `--output [-o]`    | `'dataset'`   | name of the output dataset                       |
| `--fmt [-f]`       | `'auto'`      | format of input dataset                          |
| `--ofmt [-of]`     | `'runner'`    | format of output dataset                         |
| `--emap [-em]`     | `None`        | map the elements according to a LAMMPS data file |
| `--shuffle`        | `False`       | shuffle the dataset (only for indexable formats) |
| `--seed`           | `0`           | random seed for shuffling                        |
| `--cp2k-ener`      | `'1.ener'`    | suffix for CP2K energy file                      |
| `--cp2k-cell`      | `'1.cell'`    | suffix for CP2K cell file                        |
| `--cp2k-pos`       | `'pos-1.xyz'` | suffix for CP2K position file                    |
| `--cp2k-frc`       | `'frc-1.xyz`  | suffix for CP2K force file                       |
