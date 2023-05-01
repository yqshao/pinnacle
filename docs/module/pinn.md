# PiNN workflows

## pinnTrain

The `pinnTrain` process runs a training given a dataset and a PiNN input file.
The input file can optionally be a model folder in such case the process
contines training from the last checkpoint available. For list of flags
available, see the [PiNN
documentation](https://teoroo-cmc.github.io/PiNN/master/usage/cli/train/).

### Input Channels

| Channel   | Type   | i/o            | Note                             |
|-----------|--------|----------------|----------------------------------|
| `name`    | `val`  | `in[0]`        | an id to identify the process    |
| `dataset` | `file` | `in[1]`        | a dataset recognizable by PiNN   |
| `input`   | `file` | `in[2]`        | a PiNN `.yml` input file         |
| `flag`    | `val`  | `in[3]`        | flags for `pinn train`           |
| `name`    | `val`  | `out.*[0]`     | same as input                    |
| `model`   | `file` | `out.model[1]` | Trained PiNN model               |
| `log`     | `file` | `out.log[2]`   | PiNN log (the evaluation errors) |

## pinnMD

The `pinnMD` process takes a trained model and runs a MD trajecotry, with some
limited options regarding the dynamics. The process suppors PiNN models only for
now. A list of models can be supplied to run an ensemble MD. For more complex
processes, consider writing a customized MD process and include it in your
workflow.

### Input Channels

| Channel | Type   | i/o           | Note                                                                |
|---------|--------|---------------|---------------------------------------------------------------------|
| `name`  | `val`  | `in[0]`       | an id to identify the process                                       |
| `model` | `file` | `in[1]`       | trained PiNN model, can be a list of models                         |
| `init`  | `file` | `in[2]`       | initial geometry, in any ASE recognizable format                    |
| `flags` | `val`  | `in[3]`       | a string specifying the MD simulation see [options](#options) below |
| `name`  | `val`  | `out.*[0]`    | same as input                                                       |
| `traj`  | `file` | `out.traj[1]` | output trajectory                                                   |
| `log`   | `file` | `out.log[1]`  | MD log                                                              |

### Options

The flags for pinnMD is specified in the form of `--option1 val1 --option2 val`,
the available options and default valeus are listed below.

| Option              | Default   | Note                                   |
| ------------------- | --------- | -------------------------------------- |
| `--ensemble `       | `'nvt' `  | `'npt'` or `'nvt'`                     |
| `--T `              | `273 `    | Temperature in K                       |
| `--t `              | `100 `    | Time in ps                             |
| `--dt `             | `0.5 `    | Time step in fs                        |
| `--taut `           | `100 `    | Damping factor for thermostat in steps |
| `--taup `           | `1000 `   | Damping factor for barostat in steps   |
| `--log-every `      | `5 `      | Log interval in steps                  |
| `--pressure `       | `1 `      | pressure in bar                        |
| `--compressibility` | `4.57e-5` | compressibility in bar                 |

??? "Source code"

    ```groovy
    --8<-- "nextflow/module/pinn.nf"
    ```
