# TIPS module

The `tips.nf` module contains several processes supplied by the TIPS library.

## convertDS

The `convertDS` process converts a one dataset to another. The input/output
formats are controlled by the `flags` channel.

### Channel specification

| Element     | Type   | i/o      | Note                              |
|-------------|--------|----------|-----------------------------------|
| `name`      | `val`  | `in[0]`  | an id to identify the process     |
| `input`     | `path` | `in[1]`  | input dataset                     |
| `flag`      | `val`  | `in[2]`  | flags for `tips convert`          |
| `name`      | `val`  | `out[0]` | same as input                     |
| `converted` | `path` | `out[1]` | converted dataset [`converted.*`] |

## mergeDS

The `mergeDS` process merges a number of single point calculations into one.
Note that the process also expect a `idx` element in the input channel, which
should give an index of the corresponding single point calculation, and will be
saved into the `merged.idx` file.

### Channel specification

| Element  | Type   | i/o      | Note                                         |
|----------|--------|----------|----------------------------------------------|
| `name`   | `val`  | `in[0]`  | an id to identify the process                |
| `idx`    | `val`  | `in[1]`  | indices of single point simulations          |
| `logs`   | `path` | `in[2]`  | logs from single point computations          |
| `name`   | `val`  | `out[0]` | same as input                                |
| `idx`    | `path` | `out[1]` | file that records the indices [`merged.idx`] |
| `merged` | `path` | `out[2]` | merged dataset [`merged.traj`]               |

## mixDS

The `mixDS` process takes two datasets, called `newDS` and `oldDS`, and two
flags `newFlag` and `oldFlag`, the datasets are first subsampled with
corresponding flags, and them merged together. This process is mainly used to
update a training set in an activated learning loop.

### Channel specification

| Element   | Type   | i/o      | Note                          |
|-----------|--------|----------|-------------------------------|
| `name`    | `val`  | `in[0]`  | an id to identify the process |
| `newDS`   | `path` | `in[1]`  | new dataset                   |
| `oldDS`   | `path` | `in[2]`  | old dataset                   |
| `newFlag` | `path` | `in[3]`  | subsample flag for newDS      |
| `oldFlag` | `path` | `in[4]`  | subsample flag for oldDS      |
| `name`    | `val`  | `out[0]` | same as input                 |
| `idx`     | `path` | `out[1]` | merged index (`merged.idx`)   |

## checkConverge

This workflow compares a sampled trajectories to labelled data. The output
geometry will be:

- The last frame of the trajectory if the trajectory is deemed converted;
- The first frame of the trajectory otherwise.

The convergence is controlled by the following parameters.

### Channel specification

| Element | Type   | i/o      | Note                                |
|---------|--------|----------|-------------------------------------|
| `name`  | `val`  | `in[0]`  | an id to identify the process       |
| `idx`   | `path` | `in[1]`  | index of labels in the trajectory   |
| `label` | `val`  | `in[2]`  | labelled data set                   |
| `traj`  | `val`  | `in[3]`  | sampled trajectory                  |
| `name`  | `val`  | `out[0]` | same as input                       |
| `geo`   | `path` | `out[1]` | geometry                            |
| `out`   | `val`  | `out[2]` | a string of convergence information |

### Parameters

| Parameter   | Default | Description               |
|-------------|---------|---------------------------|
| `fmaxtol`   | `2.0`   | Max error on forces       |
| `emaxtol`   | `0.02`  | Max error on energy       |
| `frmsetol ` | `0.15`  | Tolerance for force RMSE  |
| `ermsetol ` | `0.005` | Tolerance for energy RMSE |


??? "Source code"

    ```groovy
    --8<-- "nextflow/module/tips.nf"
    ```
     
