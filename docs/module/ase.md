#  ASE workflows

## aseMD

The aseMD process takes a trained model and runs a MD trajecotry, with some
limited options regarding the dynamics. For more complex processes, consider
writing a customized MD process and include it in your workflow.

### Channels

| Channel    | Format | Note                                             |
|------------|--------|--------------------------------------------------|
| (in) model | file   | a trained ANN model                              |
| (in) tag   | string | a string specifying the MD simulation see below  |
| (in) geo   | file   | initial geometry, in any ASE recognizable format |
| (out) traj | file   | output trajectory                                |

### Options

| Option    | Default | Note             |
|-----------|---------|------------------|
| -nvt/-npt | nvt     | Ensemble         |
| -T        | 273     | Temperature in K |
| -P        | 1       | Pressure in bar  |
| -dt       | 0.5     | Time step in ps  |
| -time     | 100     | Time in ps       |


??? "Source code"

    ```Python
    #!/usr/bin/env Python
    
    def __main__():
        return traj
    
    __main__()
    ```
