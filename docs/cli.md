# Command Line Interface (CLI)

*TIPS* provides a command interface as gluing utilities between different
trainer and sampler modules. 

## Commands

### convert 

Convert datasets between different formats.

#### Usage

```
$ tips convert dsfile1 [dsfile2 ...] [options]
```

#### Description

| Name, shorthand (if any) | Default      | Description                |
|--------------------------|--------------|----------------------------|
| -o                       | 'dataset.db' | name of the output dataset |
| --format, -f             | 'auto'       | format of input dataset    |
| --oformat, -of           | 'auto'       | format of output dataset   |

#### Example

From a lammps trajectory `npt.dump` to a dataset

```
$ tips convert npt.dump npt2.dump -o dataset 
```

### filter 

Convert datasets between different formats.

```
$ tips filter npt.dump
```

#### Algorithms

- **'naive'**: naive filtering given a per-structure errror or uncertaintly estimation.
- **'fps'**: furthest point sampling as outlined in [ref:]

#### Description

| Name, shorthand (if any) | Default        | Description                                    |
|--------------------------|----------------|------------------------------------------------|
| --algo, -a               | 'naive'        |                                                |
| --tol, -t                | 5%             | error tolerance (absolute value or percentage) |
| --error, -e              | 'error.dat'    | error/uncertainty label                        |
| --fingerprint, -fp       | 'fingerprints' | fingerprint tag                                |

### split

Split dataets to subsets.

#### Description

| Name, shorthand (if any) | Default          | Description                     |
|--------------------------|------------------|---------------------------------|
| --splits, -s             | 'train:8,test:2' |                                 |
| --shuffle                | true             |                                 |
| --seed                   | 0                | random seed if shufffle is used |


### reduce (to be implemented)

Dimension reduction.
