# Activated Learning

The activated learning recipe actively samples the potential. The workflow is
controlled by the several subworkflows, including the training, sampling, and
labelling processes.

## Workflow

The workflow is shown in the following diagram:

![workflow](../tikz/acle.svg)

The results of individual processes will be stored in the following directories:

- `gen${i}/emd`: (ensemble) molecular dynamics trajectories;
- `gen${i}/label`: labelled single point structures;
- `gen${i}/merge`: merged labels;
- `gen${i}/mix`: new training set for next generation;
- `gen${i}/model`: trained models.

## Parameters

The workflow supports a large number of parameters to adjust a runtime, most of
those parameters are set to reasonable defaults if you generate your workflow
from a template, see below for a detailed list for all parameters.

### Initalization Options

| Parameter      | Default                        | Description                                                         |
|----------------|--------------------------------|---------------------------------------------------------------------|
| `restart_from` | `false`                        | restart from a given generation                                     |
| `restart_conv` | `false`                        | restart from a converged model (model will be retrained if `false`) |
| `proj `        | `'acle'`                       | folder for storing results                                          |
| `init_geo `    | `'input/geo/*.xyz'`            | inital geometries for sampling                                      |
| `init_model `  | `'input/pinn/pinet-adam.yml'`  | initial model or model parameters                                   |
| `init_ds `     | `'input/ds/init-ds.{yml,tfr}'` | initial dataset                                                     |
| `init_time `   | `1.0`                          | sampling time scale in ps                                           |
| `init_steps `  | `100000`                       | training steps for initial model                                    |
| `ens_size `    | `1`                            | size of the model ensemble (for ensemble MD sampling)               |
| `geo_size `    | `6`                            | size of the starting geometries for sampling                        |

### Iteration Options

| Parameter      | Default            | Description                                         |
|----------------|--------------------|-----------------------------------------------------|
| `retrain_step` | `50000 `           | number of retrain steps per generation              |
| `label_flags ` | `'--nsample 50' `  | selection rule for the data to label                |
| `old_flag `    | `'--nsample 2700'` | selection rule for the old dataset                  |
| `new_flag `    | `'--nsample 300'`  | selection rule for the new dataset                  |
| `sp_points `   | `50`               | number of single points for each sampled trajectory |
| `acc_fac `     | `2.0 `             | factor to acceralate/slow down the sampling         |
| `min_time `    | `1.0 `             | minimal timescale for sampling                      |
| `emaxtol `     | `0.020 `           | toleranace for max error error                      |
| `ermsetol `    | `0.005 `           | toleranace for energy RMSE                          |
| `fmaxtol `     | `0.800 `           | toleranace for max force component error            |
| `frmsetol `    | `0.200 `           | toleranace for force RMSE                           |

??? "Source code"

    ```groovy
    --8<-- "nextflow/acle.nf"
    ```
