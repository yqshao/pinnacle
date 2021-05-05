#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.publishDir      = 'tips'
params.publishMode     = 'link'
include {fileList; getParams} from "$moduleDir/utils"

defaults = [:]
defaults.subDir       = '.'
defaults.ds           = null
defaults.inp          = ''
defaults = getParams(defaults, params)
process tipsFilter {
    publishDir {"$params.publishDir/$setup.subDir"}, mode: params.publishMode
    label 'tips'

    input:
    tuple val(meta), val(inputs)

    output:
    tuple val(meta), path('output.xyz')

    script:
    setup = getParams(defaults, inputs)
    """
    #!/usr/bin/env bash
    tips filter ${fileList(setup.ds)} $setup.inp -o output -of xyz
    """

    stub:
    setup = getParams(defaults, inputs)
    """
    #!/usr/bin/env bash
    touch output.xyz
    """
}
