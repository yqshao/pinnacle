#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.publishDir      = 'pinn'
params.publishMode     = 'link'

include {fileList; getParams} from "$moduleDir/utils"

trainDflts = [:]
trainDflts.subDir            = '.'
trainDflts.inp               = null
trainDflts.ds                = null
trainDflts.seed              = 0
trainDflts.maxSteps          = '1000000'
trainDflts.genDress          = true
trainDflts.pinnCache         = 'True'
trainDflts.pinnBatch         = '10'
trainDflts.pinnCkpts         = '1'
trainDflts.pinnShuffle       = '500'
trainDflts.pinnLogSteps      = '1000'
trainDflts.pinnCkptSteps     = '10000'
trainDflts = getParams(trainDflts, params)
process pinnTrain {
    label 'pinn'
    publishDir {"$params.publishDir/$setup.subDir"}, mode: params.publishMode

    input:
    tuple val(meta), val(inputs)

    output:
    tuple val(meta), path('model/', type:'dir')

    script:
    setup = getParams(trainDflts, inputs)
    """
    tips convert ${fileList(setup.ds)} -o 'train:8,eval:2' -of pinn --seed $setup.seed
    if [ ! -f ${file(setup.inp)}/params.yml ];  then
        mkdir -p model; cp ${file(setup.inp)} model/params.yml
    else
        cp -r ${file(setup.inp)} model; rm -r model/events* model/eval
    fi
    pinn_train --model-dir='model' --params-file='model/params.yml'\
        --train-data='train.yml' --eval-data='eval.yml'\
        --cache-data=$setup.pinnCache\
        --batch-size=$setup.pinnBatch\
        --shuffle-buffer=$setup.pinnShuffle\
        --train-steps=$setup.maxSteps\
        --max-ckpts=$setup.pinnCkpts\
        --log-steps=$setup.pinnLogSteps\
        --ckpt-steps=$setup.pinnCkptSteps\
        ${setup.genDress? "--regen-dress": ""}
    """

    stub:
    setup = getParams(trainDflts, inputs)
    """
    mkdir model
    """
}

labelDflts = [:]
labelDflts.subDir       = '.'
labelDflts.inp          = null
labelDflts.ds           = null
labelDflts = getParams(labelDflts, params)
process pinnLabel {
    label 'pinn'
    publishDir {"$params.publishDir/$setup.subDir"}, mode: params.publishMode

    input:
    tuple val(meta), val(inputs)

    output:
    tuple val(meta), path('output.xyz')

    script:
    setup = getParams(labelDflts, inputs)
    """
    #!/usr/bin/env python3
    import pinn, yaml, os
    import tensorflow as tf
    from ase import units
    from ase.io import read, write
    os.symlink("${file(setup.inp)}", "model")
    calc = pinn.get_calc("model")
    traj = read("${file(setup.ds)}", index=':')
    with open('output.xyz', 'w') as f:
        for atoms in traj:
            atoms.wrap()
            atoms.set_calculator(calc)
            atoms.get_potential_energy()
            write(f, atoms, format='extxyz', append='True')
    """

    stub:
    setup = getParams(labelDflts, inputs)
    """
    #!/usr/bin/env bash
    touch output.xyz
    """
}

sampleDflts = [:]
sampleDflts.subDir       = '.'
sampleDflts.inp          = null
sampleDflts.init         = null
sampleDflts.pinnEnsemble = 'NPT'
sampleDflts.pinnCopies   = 1
sampleDflts.pinnStep     = 0.5  // in fs
sampleDflts.pinnTime     = 5    // in ps
sampleDflts.pinnEvery    = 0.01 // in ps
sampleDflts.pinnTaut     = 100  // in steps
sampleDflts.pinnTaup     = 1000 // in steps
sampleDflts.pinnTemp     = 330  // in K
sampleDflts.pinnPress    = 1    // in bar
sampleDflts.pinnCompress = '4.57e-5' // in bar^-1
sampleDflts = getParams(sampleDflts, params)
process pinnSample {
    label 'pinn'
    publishDir {"$params.publishDir/$setup.subDir"}, mode: params.publishMode

    input:
    tuple val(meta), val(inputs)

    output:
    tuple val(meta), path('output.xyz')

    script:
    setup = getParams(sampleDflts, inputs)
    """
    #!/usr/bin/env python3
    import pinn, os
    import numpy as np
    import tensorflow as tf
    from ase import units
    from ase.io import read, write
    from ase.io.trajectory import Trajectory
    from ase.md import MDLogger
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md.nptberendsen import NPTBerendsen
    from ase.md.nvtberendsen import NVTBerendsen

    os.symlink("${file(setup.inp)}", "model")
    calc = pinn.get_calc("model")
    ensemble = "$setup.pinnEnsemble"
    for seed in range($setup.pinnCopies):
        rng = np.random.default_rng(seed)
        atoms = read("${file(setup.init)}")
        atoms.set_calculator(calc)
        MaxwellBoltzmannDistribution(atoms, $setup.pinnTemp*units.kB, rng=rng)
        dt = 0.5 * units.fs
        steps = int($setup.pinnTime*1e3*units.fs/dt)
        if ensemble=='NPT':
            dyn = NPTBerendsen(atoms, timestep=dt, temperature=$setup.pinnTemp, pressure=$setup.pinnPress,
                               taut=dt*$setup.pinnTaut, taup=dt*$setup.pinnTaup, compressibility=$setup.pinnCompress)
        elif ensemble=='NVT':
            dyn = NVTBerendsen(atoms, timestep=dt, temperature=$setup.pinnTemp, taut=dt * 100)
        else:
            raise NotImplementedError(f"Unkown ensemble {ensemble}")
        interval = int($setup.pinnEvery*1e3*units.fs/dt)
        dyn.attach(MDLogger(dyn, atoms, 'output.log', mode="a"), interval=interval)
        dyn.attach(Trajectory('output.traj', 'a', atoms).write, interval=interval)
        try:
            dyn.run(steps)
        except:
            pass
    traj = read('output.traj', index=':')
    [atoms.wrap() for atoms in traj]
    write('output.xyz', traj)
    """

    stub:
    setup = getParams(sampleDflts, inputs)
    """
    #!/usr/bin/env bash
    touch output.xyz
    """
}

workflow {
    pinnTrain(Channel.of([:]))
}
