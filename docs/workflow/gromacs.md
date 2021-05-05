# tips/gromacs

This module implements the `sampler` and `labeller` process with the
[Gromacs](https://www.gromacs.org/) molecular dynamics code. All Gromacs
processes shares a `gmxTop` input that specifies the topology file required to
generate the `.mdp` input. With in the process, the structure or trajectories
are converted to the `.pdb` format and then used to generate the `.tpr` files.

## gromacsSample

| Inputs | Default | Description |
|--------|---------|-------------------|
| inp    | `null`  | Gromacs mdp file  |
| init   | `null`  | initial structure |
| gmxTop | `null`  | Gromacs top file |

Outputs the sampled trajectory.

## gromacsLabel
Label a dataset with Gromacs

| Inputs | Default | Description      |
|--------|---------|------------------|
| inp    | `null`  | Gromacs mdp file |
| ds     | `null`  | dataset file     |
| gmxTop | `null`  | Gromacs top file |

Outputs the labelled dataset.
