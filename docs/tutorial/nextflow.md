# More on nextflow

## The deisgn of PiNNAcLe workflows

The workflows in PiNNAcLe follows the
[DSL2](https://www.nextflow.io/docs/latest/dsl2.html#) syntax of Nextflow, which
allows processes and sub-workflows to be reused and imported as modules. In
additon, the basic processes (trainers, samplers, and labellers) shares a
similar structures, so that more complex workflows can be adapted to different
implementations easily, taking `pinnMD` sampler as an example

```groovy
process pinnMD {
  tag "$name" // (1)
  label 'pinn' 
  publishDir "$params.publish/$name" // (2)

  input:
    tuple val(name), path(model,stageAs:'model*'), path(init), val(flags) // (3)

  output:
    tuple val(name), path('pinnmd.traj'), emit: traj // (4)
    tuple val(name), path('pinnmd.log'), emit: log

    ....
}
```

1. All PiNNAcLe workflows accepts tuples as input, where the "name" is a unique
   identifier of a process containning necessary metadata.
2. PiNNAcLe workflows typically uses the `publish` parameter to specify the "parent"
   folder of the same process, while different instantiation are placed into
   subfolders named after "name".
3. Processes of the same "class" shares a same input structure, here, a MD
   sampler always takes a model and an initial structure as input. The
   implementation specific settings will be merged into the input `flags`.
4. Processes and have multiple outputs, again, similar processes will have the
   same naming pattern for outputs.

## Patterns in complex workflows

### Changing publish folder

When implementing a complex workflow it often make sense to arrange the
different results with into individual subfolders, the following snippet
pubishes the results into subfolders of the `publish` parameter as specified in
the main script.

```groovy
include { cp2k } from './tips/cp2k.nf' addParams(publish: "$params.proj/cp2ksp")
include { pinnTrain } from './tips/pinn.nf' addParams(publish: "$params.proj/models")
```

### Branching and merging

In the order to parallelize the processes and combine them latter, PiNNAcLe
always use the "name" as an identifier, the branched processes will be labelled
as `name/branchname`, and they might be combined later with the
[groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple)
operator. (we use regular expressions to match naming patterns, when unsure,
check [regex101](https://regex101.com))

```groovy
ch_models \
  | combine (ch_init) \
  | map {gen, models, init -> \
         ["model.baseName/$init.baseName", model, init, params.md_flags]} \
  | pinnMD // (1) 

pinnMD.out.traj \
  | map {name, traj -> \ 
         [(name=~/(.+)\/.+/)[0][1], traj]} \
  | groupTuple  // (2)
```

1. the above pipeline iterates over combinations of models and initial geometries.
2. the lower part collects the trajectories from the same model.

