# tips/pinn

This module implements the `trainer`, `sampler` and `labeller` processes with
the [PiNN](https://github.com/Teoroo-CMC/pinn) package. The `inp` should be a
PiNN model directory, for the `pinnTrain` process, the inp can 

## pinnTrain

| Inputs        | Default   | Description                                |
|---------------|-----------|--------------------------------------------|
| inp           | `null`    | PiNN model or params file                  |
| ds            | `null`    | trianing set                               |
| seed          | `0`       | random seed for generating training set    |
| maxSteps      | `1000000` | max iteration [in steps]                   |
| genDress      | `true`    | regenerate the atomic dress                |
| pinnCache     | `'True'`  | cache preprocessed dataset during training |
| pinnBatch     | `10`      | batch size                                 |
| pinnCkpts     | `1`       | max number of checkpoints to save          |
| pinnShuffle   | `500`     | shuffle buffer size                        |
| pinnLogSteps  | `1000`    | save summary every x steps                 |
| pinnCkptSteps | `10000`   | save checkpoint every x steps              |

Outputs the final model folder with checkpoints.

## pinnLabel
Label a dataset with PiNN

| Inputs | Default | Description  |
|--------|---------|--------------|
| inp    | `null`  | PiNN model   |
| ds     | `null`  | trianing set |

Outputs the labelled dataset.

## pinnSample

| Inputs       | Default   | Description                            |
|--------------|-----------|----------------------------------------|
| inp          | `null`    | PiNN model                             |
| init         | `null`    | initial Structure                      |
| pinnEnsemble | `'NPT'`   | `'NPT'` or `NVT`                       |
| pinnCopies   | `1`       | resample the traj with different seeds |
| pinnStep     | `0.5`     | timestep for MD simulation             |
| pinnTime     | `5`       | tot sample time [ps]                   |
| pinnEvery    | `0.01`    | sample every [ps]                      |
| pinnTaut     | `100`     | sample every [steps]                   |
| pinnTaup     | `1000`    | sample every [steps]                   |
| pinnTemp     | `330`     | sample at temperature [K]              |
| pinnPress    | `1`       | sample at pressure [bar]               |
| pinnCompress | `4.578-5` | compressibility for barastat [bar^-1]  |

Outputs the sampled dataset.
