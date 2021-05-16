# tips/gromacs

This module implements the `sampler` and `labeller` process with the
[Gromacs](https://www.gromacs.org/) molecular dynamics code. All Gromacs
processes shares a `gmxTop` input that specifies the topology file required to
generate the `.tpr` input. The `mdp` files used in both process can be similar,
the only difference between them is that `gromacsLabel` uses the `--rerun`
option of `gmx mdrun` (you might want to set to set the `nst` options to `1` in
the input for labelling process).

In both processes, the structure or trajectories are converted to the `.pdb`
format and then used to generate the `.tpr` files. The same `.pdb` file is used
as the `top` file when converting the outputs to datasets. The files are
converted using the MDAnalysis code and the energy units will be `kcal/mol`.

## gromacsSample

| Inputs | Default | Description       |
|--------|---------|-------------------|
| inp    | `null`  | Gromacs mdp file  |
| init   | `null`  | initial structure |
| gmxTop | `null`  | Gromacs top file  |

Outputs the sampled trajectory.

## gromacsLabel
Label a dataset with Gromacs

| Inputs | Default | Description      |
|--------|---------|------------------|
| inp    | `null`  | Gromacs mdp file |
| ds     | `null`  | dataset file     |
| gmxTop | `null`  | Gromacs top file |

Outputs the labelled dataset.
