#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
params.publishDir      = 'gromacs'
params.publishMode     = 'link'
include {getParams} from "$moduleDir/utils"

labelDflts = [:]
labelDflts.subDir       = '.'
labelDflts.inp          = null
labelDflts.ds           = null
labelDflts.gmxInp       = null
labelDflts = getParams(labelDflts, params)
process labeller {
    label 'gromacs'
    publishDir {"$params.publishDir/$setup.subDir"}, mode: params.publishMode

    input:
    tuple val(meta), val(inputs)

    output:
    tuple val(meta), path('output.xyz')

    script:
    setup = getParams(labelDflts, inputs)
    """
    #!/usr/bin/env bash
    ln -s ${file(setup.gmxTop)} input.init
    ln -s ${file(setup.inp)} input.mdp
    tips convert ${file(setup.ds)} -o input -of pdb
    gmx grompp -f input.mdp -c input.pdb -p input.top -o input.tpr
    gmx mdrun -s input.tpr -rerun input.pdb --deffnm output
    gmx traj-conv -f output.trr -o output.pdb
    tips convert output.pdb output.xyz
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
sampleDflts = getParams(sampleDflts, params)
process sampler {
    label 'gromacs'
    publishDir {"$params.publishDir/$setup.subDir"}, mode: params.publishMode

    input:
    tuple val(meta), val(inputs)

    output:
    tuple val(meta), path('output.xyz')

    script:
    setup = getParams(sampleDflts, inputs)
    """
    #!/usr/bin/env bash
    tips convert ${file(setup.init)} -o input -of pdb
    ln -s ${file(setup.gmxTop)} input.init
    ln -s ${file(setup.inp)} input.mdp
    gmx grompp -f input.mdp -c input.pdb -p input.top -o input.tpr
    gmx mdrun -s input.tpr --deffnm output
    gmx traj-conv -f output.trr -o output.pdb
    tips convert output.pdb output.xyz
    """

    stub:
    setup = getParams(sampleDflts, inputs)
    """
    #!/usr/bin/env bash
    touch output.xyz
    """
}
