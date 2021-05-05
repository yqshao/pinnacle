#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.publishDir      = 'lammps'
params.publishMode     = 'link'
include {getParams} from "$moduleDir/utils"

labelDflts = [:]
labelDflts.subDir       = '.'
labelDflts.inp          = null
labelDflts.ds           = null
labelDflts.lmpInit      = null
labelDflts.lmpData      = null
labelDflts.lmpSetting   = null
labelDflts.lmpEmap      = ''
labelDflts.lmpUnits     = 'real'
labelDflts.lmpCmd       = 'lmp_mpi'
labelDflts = getParams(labelDflts, params)
process lammpsLabel {
    label 'lammps'
    publishDir {"$params.publishDir/$setup.subDir"}, mode: params.publishMode

    input:
    tuple val(meta), val(inputs)

    output:
    tuple val(meta), path('output.xyz')

    script:
    setup = getParams(labelDflts, inputs)
    remap = setup.lmpEmap.tokenize(',').collect{it.tokenize(':').reverse().join(':')}.join(',')
    """
    #!/usr/bin/env bash
    ln -s ${file(setup.lmpInit)} input.init
    ln -s ${file(setup.lmpData)} input.data
    ln -s ${file(setup.lmpSetting)} input.setting
    tips convert ${file(setup.ds)} --emap '$remap' -o input -of lammps
    $setup.lmpCmd -in ${file(setup.inp)} || echo LAMMPS aborted
    sed -i '/WARNING/d' output.log
    tips convert output.dump --log output.log -o output -of xyz\
        --units $setup.lmpUnits --emap '$setup.lmpEmap'
    """

    stub:
    setup = getParams(labelDflts, inputs)
    """
    #!/usr/bin/env bash
    touch output.xyz
    """
}

sampleDflts = [:]
sampleDflts.subDir      = '.'
sampleDflts.inp         = null
sampleDflts.init        = null
sampleDflts.lmpInit     = null
sampleDflts.lmpData     = null
sampleDflts.lmpSetting  = null
sampleDflts.lmpUnits    = 'real'
sampleDflts.lmpEmap     = ''
sampleDflts.lmpCmd      = 'lmp_mpi'
sampleDflts = getParams(sampleDflts, params)
process lammpsSample {
    label 'lammps'
    publishDir {"$params.publishDir/$setup.subDir"}, mode: params.publishMode

    input:
    tuple val(meta), val(inputs)

    output:
    tuple val(meta), path('output.xyz')

    script:
    setup = getParams(sampleDflts, inputs)
    """
    #!/usr/bin/env bash
    ln -s ${file(setup.lmpInit)} input.init
    ln -s ${file(setup.lmpData)} input.data
    ln -s ${file(setup.lmpSetting)} input.setting
    $setup.lmpCmd -in ${file(setup.inp)} || echo LAMMPS aborted
    sed -i '/WARNING/d' output.log
    tips convert output.dump --log output.log \
      --units $setup.lmpUnits --emap '$setup.lmpEmap' -o output -of xyz
    """

    stub:
    setup = getParams(sampleDflts, inputs)
    """
    #!/usr/bin/env bash
    touch output.xyz
    """
}
