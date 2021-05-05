# Implementation patterns

The TIPS workflows are implemented in [Nextflow](https://www.nextflow.io/), if
you are interested in using TIPS, it would be good if you are familiar with
basic concepts in Nextflow such are data channels, `params`, and basic options 
like running and importing workflows.

This document explains the patterns that TIPS follows in implementing AML
workflows. A typical process in TIPS looks like this:

```groovy linenums="1"
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
params.publishDir   = 'myoutput'
params.publishMode  = 'link'
include {getParams} from './tips/utils'
defaults            = [:]
defaults.inp        = null
defaults.ds         = null
defaults.myParam    = 'myParam'
defaults = getParams(defaults, params)
process myworkflow {
    publishDir {"$params.publishDir"}, mode: params.publishMode
    input: tuple val(meta), val(input)
    output: tuple val{meta}, file("output.xyz") 
    script:
    setup = getParams(defaults, inputs)
    """
    my_code $setup.ds $setup.input $setup.myParam > output.xyz
    """
}

workflow {labeller([null, [:]])}
```

## Input/output

As mentioned in the [overview](overview.md#inclusion), all workflows in TIPS has
two inputs, a meta data and a "real" input. When the workflow is included as a
sub-workflow, the metadata will be set to identify the result (e.g. in
`tips/explore` this includes the iteration and random seed for each model.)

More specifically, the `inputs` should be a Groovy map (dictionary) that
contains all the options required to run a process. The reason of using a single
map for all inputs it to allow the same process to be used with different types 
of input channels.

## Default parameters 

The defaults variable holds all the variables that may be updated by `inputs`,
and their default values.You might notice that in L11 and L16, the
`getParams(source, updates)` are used twice to update the defaults, first with
the `params` and then with the `inputs`. This means the inputs to your process
can be specified with CLI options and the data channel. After the updates, the
`setup` variable should be used in the script.


## Running and caching

The example module provides a default workflow that invokes the process with an
empty map. This allows you to run the workflow as:

```
$ nf run mymodule.nf --inp inp_file --ds ds_file
```

Alternatively, within a wokflow:
```groovy
include {myworkflow} from mymodule
myworkflow([inp: 'input_file', ds: 'ds_file'])
```

The difference between the two ways is that in the later case the variables will
be considered as inputs by Nextflow, and used to compute the [cache
keys](https://www.nextflow.io/docs/latest/process.html#cache).

This should not cause a big difference if your input in used in the script.
However, if you have some variable that does not change the outcome of the
workflow (e.g. the publish options here), you might want to exclude it from you
defaults map and only put it as a key of `params` (so that you can use the
`-resume` option after changing `publishDir`).
