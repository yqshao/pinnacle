# Custom process

Some basics about processes in TIPS is introduced in [patterns](patterns.md). In
this document there are some more suggestions if you are to develop a custom
process that fits into the TIPS suite.

## subDir

`subDir` is a special input key that TIPS uses in strategies to control the
output folder in the sub-workflows. For example, the `explore` strategy will
organize the training and sampling results according to iterations models. In
you custom module, you might want to allow this behavior by setting the
publishDir to `$params.publishDir/$setup.subDir`.

## File formats

At this point most processes outputs are converted to `.xyz` formats, which is
quite readable handy for debugging. Since in principle TIPS will be able to
handle different formats of outputs, consider converting them to a specific
format before running your code, and converting them back to `.xyz` after you
have done.

## stub

As of 20.11.0-edge, Nextflow supports a [`stub`
directive](https://www.nextflow.io/docs/edge/process.html#stub) which allows for
providing a dummy command for each process. All TIPS processes comes with this
directive, if you are working on a workflow based on TIPS, you can test it with
a dry run `nf run ... -stub`.

## Template
Below is a template for implementing a custom labeller for TIPS

```groovy
params.publishDir      = 'mymodule'
params.publishMode     = 'link'

include {getParams} from "$tipdDir/utils"

defaults = [:]
defaults.subDir  = '.'
defaults.inp     = null
defaults.ds      = null
defaults.myParam = null
defaults = getParams(defaults, params)
process labeller {
    publishDir {"$params.publishDir/$setup.subDir"}, mode: params.publishMode

    input:
    tuple val(meta), val(input)

    output:
    tuple val{meta}, file("output.xyz") 

    script:
    setup = getParams(defaults, inputs)
    """
    my_code $setup.ds $setup.input $setup.myParam > output.xyz
    """
    
    stub:
    setup = getParams(defaults, inputs)
    """
    touch output.xyz
    """
}

workflow {
    labeller([null, [:]])
}
```

