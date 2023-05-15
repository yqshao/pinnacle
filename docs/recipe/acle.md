# Activated Learning

The activated learning recipe actively samples a potential energy surface. The
workflow is controlled by several subworkflows, including the training,
sampling, and labelling processes. The overall workflow is shown in the
following diagram:

![workflow](../tikz/acle.svg)

The workflow can be used as either an entrypoint or a subworkflow. Some
parameters that set up the initial datasets and models are taken by the
entrypoint only, see below for usage and tables of parameters.

## Entrypoint 

??? Code "Convert initial dataset and geometry from a ASE trajectory"

    As a example, consider a trajectory file readable by ASE as input, 
    the `tips convert` CLI tool can be used get the initial dataset and
    geometries:
    
    ```bash
    # generating my_ds.{yml,tfr} 
    tips convert -f asetraj data.traj -of pinn my_ds
    # generating the initial geometry
    tips convert -f asetraj data.traj --subsample uniform -of idx.xyz
    ```

To run the workflow as an entrypoint (**single quotes are necessary**):

```bash
nextflow run main.nf -entry acle --init_ds 'myds.{yml,tfr}' --init_geo '*.xyz' ...
```

It is possible to restart a project from a ceratin generation, while keeping the
folder structure:

```bash
nextflow run main.nf -entry acle --restart_from 30 --restart_conv true
```

When restarted in the above way the `init_*` parameters will be ignored. This
method is mainly for small changes of the sampled ensemble, e.g., an *ad hook*
change of temperature. In cases where this is not enough, it is advisable to
rerun the workflow under a different `--proj`.

### Parameters

| Parameter      | Description                                                         | Default                      |
|----------------|---------------------------------------------------------------------|------------------------------|
| `init_geo`     | inital geometries for sampling                                      | `input/geo/*.xyz`            |
| `init_model`   | initial model or model parameters                                   | `input/pinn/pinet-adam.yml`  |
| `init_ds`      | initial dataset                                                     | `input/ds/init-ds.{yml,tfr}` |
| `init_time`    | sampling time scale in ps                                           | `1.0`                        |
| `init_steps`   | training steps for initial model                                    | `100000`                     |
| `restart_from` | restart from a given generation                                     | `false`                      |
| `restart_conv` | restart from a converged model (model will be retrained if `false`) | `false`                      |

## Subworkflow

### Input/Output Channels


| Channel    | I/O[idx] | Type   | Description                                 |
|------------|----------|--------|---------------------------------------------|
| `gen`      | `in[0]`  | `val`  | generation of the model                     |
| `geo`      | `in[1]`  | `file` | initial geometry                            |
| `ds`       | `in[2]`  | `file` | training dataset                            |
| `steps`    | `in[3]`  | `val`  | training steps                              |
| `time`     | `in[4]`  | `val`  | sampling timescale                          |
| `converge` | `in[5]`  | `val`  | whether the input model is deemed converged |

!!! Code ""

    The AcLe subworkflow is a recursive workflow, and the input and output
    shares the same data structure.

### Parameters

| Parameter       | Description                                         | Default                                                                                |
|-----------------|-----------------------------------------------------|----------------------------------------------------------------------------------------|
| `proj`          | folder for storing results                          | `acle`                                                                                 |
| `ref`           | reference calculation module                        | `dftb`                                                                                 |
| `ref_inp`       | reference input file                                | `input/dftb/xtb.py`                                                                    |
| `mpl`           | machine learning potential module                   | `pinn`                                                                                 |
| `train_flags`   | mlp training flags                                  | `--log-every 10000 --ckpt-every 100000 --batch 1 --max-ckpts 1 --shuffle buffter 3000` |
| `train_init`    | mlp training flags                                  | `--init`                                                                               |
| `max_gen`       | maximal number of generations                       | `40`                                                                                   |
| `min_time`      | minimal timescale for sampling                      | `1.0`                                                                                  |
| `max_time`      | maximal timescale for sampling                      | `1000.0`                                                                               |
| `md_flags`      | flags for md sampling, see ase module for details   | `--ensemble nvt --dt 0.5 --log-every 100 --T 340`                                      |
| `collect_flags` | collection flags for the data to label              | `-f asetraj --subsample uniform --nsample 10 -of idx.xyz -o ds`                        |
| `sp_points`     | number of single points per sampled trajectory      | `10`                                                                                   |
| `old_flag`      | selection rule for the old dataset                  | `--nsample 2700`                                                                       |
| `new_flag`      | selection rule for the new dataset                  | `--nsample 300`                                                                        |
| `sp_points`     | number of single points for each sampled trajectory | `50`                                                                                   |
| `emaxtol`       | toleranace for max error error                      | `0.020`                                                                                |
| `ermsetol`      | toleranace for energy RMSE                          | `0.005`                                                                                |
| `fmaxtol`       | toleranace for max force (component) error          | `0.800`                                                                                |
| `frmsetol`      | toleranace for force (component) RMSE               | `0.200`                                                                                |
| `retrain_step`  | number of retrain steps per generation              | `100000`                                                                               |
| `acc_fac`       | factor to acceralate the sampling                   | `2.0`                                                                                  |
| `brake_fac`     | factor to slow down the sampling                    | `1.0`                                                                                  |

---

??? "Source Code: nextflow/acle.nf"

    ```groovy
    --8<-- "nextflow/acle.nf"
    ```
