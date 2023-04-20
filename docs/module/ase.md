# ASE workflows

## aseMD

The aseMD process takes a trained model and runs a MD trajecotry, with some
limited options regarding the dynamics. The process suppors PiNN models only for
now. A list of models can be supplied to run an ensemble MD. For more complex
processes, consider writing a customized MD process and include it in your
workflow.

### Channels

| Channel    | Type | Note                                                                |
| ---------- | ---- | ------------------------------------------------------------------- |
| (in) name  | val  | an id to identify the process                                       |
| (in) model | file | trained ANN model, can be a list of models                          |
| (in) flags | val  | a string specifying the MD simulation see [options](#options) below |
| (in) init  | file | initial geometry, in any ASE recognizable format                    |
| (in) pub   | val  | path to publish the output trajectory                               |
| (out) traj | file | output trajectory                                                   |
| (out) log  | file | MD log                                                              |

### Options

The MD specifics for aseMD is specified as a string of flags in the form of
`--option1 val1 --option2 val`, the available options and default valeus are down below.

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
    --8<-- "nextflow/module/ase.nf"
    ```
