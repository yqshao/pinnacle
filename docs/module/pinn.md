# PiNN workflows

## pinnTrain

The `pinnTrain` process runs a training given a dataset and a PiNN input file.
The input file can optionally be a model folder in such case the process
contines training from the last checkpoint available. For list of flags
available, see the [PiNN
documentation](https://teoroo-cmc.github.io/PiNN/master/usage/cli/train/).

### Channels

| Channel      | Type | Note                                     |
| ------------ | ---- | ---------------------------------------- |
| (in) name    | val  | an id to identify the process            |
| (in) dataset | file | a dataset recognizable by PiNN           |
| (in) input   | file | a PiNN `.yml` input file                 |
| (in) flag    | val  | flags for `pinn train`                   |
| (out) model  | file | Trained PiNN model                       |
| (out) log    | file | PiNN log (the evaluation errors)         |

??? "Source code"

    ```groovy
    --8<-- "nextflow/module/pinn.nf"
    ```
