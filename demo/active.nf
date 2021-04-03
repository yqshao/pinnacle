#!/usr/bin/env nextflow

// Demo active learning workflow with tips (processes indicated as boxes)
// The workflow is analogous to that of:
// C. Schran, J. Behler & D. Marx, J. Chem. Theory Comput., 2019, 16, 88-99.
//
//
// -------------------- Phase 2 Exploration by Qbc -----------------------------
//                                   ┌────────┐
//                            ┌─────►│ml_label├─►ds2Qbc
//                            │      └────────┘    │
//                            │           ▲        ▼
//             ┌───────┐      │      ┌────┴────┐ ┌───┐y
//     initDs─►│trainer├─► mlModel──►│ml_sample│ │qbc├─►finalModel
//             └───────┘             └─────────┘ └─┬─┘
//                 ▲                               │n
//                 │      ┌──────────┐             │
//                augDs◄──┤base_label│◄─ds2Label◄──┘
//                        └──────────┘
// -----------------------------------------------------------------------------

// Parameters (fixed or from input)
// Lammps Setup (unit: real):
// - InitDs: 100ps, dump every 0.1 ps (1000 points)
// - FF: water (SPC/Fw) + NaCl (Joung-Cheathem)
// - Ensemble: NVT (CSVR barostat)
// ASE Sampling:
// - TimeStep: 0.5 fs, dump every step
// - Ensemble: NPT at 1bar/330K (berendsen barostat)
// PiNN Training:
// - PiNet Parameters taken from the PiNN paper:
//   Y. Shao, M. Hellström, P. D. Mitev, L. Knijff & C. Zhang,
//   J. Chem. Inf. Model., 2020, 60, 3.

params.baseDir =  './'       //base directory for outputs
// model related
params.maxIter = 30          //no. generations
params.modelSeeds = 2        //no. models
params.modelParams = 'inputs/pinet.yml'
params.trainSteps = 200000    //steps for initDs
params.retrainSteps = 100000  //steps for each augDs
// sampling settings
params.init = 'inputs/init.{lmp,geo}'
params.labeller = 'inputs/{label.lmp,init.geo}'
params.sampleInit = 'inputs/init.xyz'
params.sampleTime = 0.1      //resample time [ps]
params.sampleInterv = 0.005  //resample every [ps]
params.sampleSeeds = 10      //resample seeds
// filters
params.qbcFilter = ''          // do not filter here
params.labelFilter = '-vmax "e:-20" -amax "f:10"'
params.restartFilter = '-vmax "e:-42" -amax "f:4"'

baseDir = params.baseDir
// Create the output channels
initDs = Channel.create()   // for training sets (iter, ds)
trainDs = Channel.create()   // for training sets (iter, ds)
models = Channel.create()   // for trained models (iter, seed, modelDir)
trajs = Channel.create()    // for trajs
ckpts = Channel.create()
restart = Channel.create()  // restart for sampling
ds4label = Channel.create() // ds to be labelled
aug4combine = Channel.create() // ds to be labelled

// Connecting the channels
Channel.of(1..params.modelSeeds)
    .map{[1, it, file(params.modelParams)]}
    .mix(ckpts)
    .set{ckpt4train}
Channel.of([2, file(params.sampleInit)])
    .mix(restart.map{iter,restart->[iter+1,restart]})
    .until{it[0]>params.maxIter}
    .set{init4sample}
initDs.map{[1, it]}
    .mix(trainDs)
    .tap{ds4combine}
    .combine(Channel.of(1..params.modelSeeds)).map{it -> it[[0,2,1]]}
    .tap{ds4train}
models
    .map{iter, seed, model -> [iter+1, seed, model]}
    .until{it[0]>params.maxIter}
    .tap(ckpts)
    .tap{finalModel}
    .branch {
        md: it[1]==1
            return [it[0], it[2]]
        other: true}
    .set{model4qbc}
trajs
    .tap{qbcRef}
    .combine(Channel.of(2..params.modelSeeds)).map{it -> it[[0,2,1]]}
    .join(model4qbc.other, by: [0,1])
    .set{qbcInp}
ds4combine
    .map{iter,ds -> [iter+1, ds]}
    .until{it[0]>params.maxIter}
    .join(aug4combine)
    .map{iter, old, aug -> [iter, ([]<<old<<aug).flatten()]}
    .tap(trainDs)

// Inputs for processes
trainInp = ckpt4train.join(ds4train, by: [0,1])
sampleInp = init4sample.join(model4qbc.md)

// Proceses
process kickoff {
    publishDir "$baseDir/datasets/", pattern: '{train,test}_1.xyz', mode: 'link'
    publishDir "$baseDir/trajs/iter1", pattern: 'label.xyz', mode: 'link'
    label 'lammps'

    input:
    path inp from Channel.fromPath(params.init).collect()

    output:
    path "train_1.xyz"  into initDs
    path "test_1.xyz"
    path "label.xyz"

    """
    mpirun -np $task.cpus lmp_mpi -in init.lmp
    tips convert prod.dump --log prod.log --units real\
         --emap '1:1,2:8,3:11,4:17' -o label -of xyz
    tips split label.xyz -s 'train_1:9,test_1:1' -of xyz
    """
}

process trainner {
    publishDir "$baseDir/models/iter$iter/seed$seed", mode: 'link'
    stageInMode 'copy'
    label 'pinn'

    input:
    tuple val(iter), val(seed), path(model), path(dataset, stageAs: 'ds/*') from trainInp

    output:
    tuple val(iter), val(seed), path('model/', type:'dir') into models

    script:
    """
    tips split ds/*.xyz  -s 'train:8,eval:2' --seed $seed
    [ ! -f model/params.yml ] && { mkdir -p model; cp $model model/params.yml; }
    pinn_train --model-dir='model' --params-file='model/params.yml'\
        --train-data='train.yml' --eval-data='eval.yml'\
        --cache-data=True --batch-size=1 --shuffle-buffer=500\
        --train-steps=${params.trainSteps+(iter-1)*params.retrainSteps}\
        ${iter==1? "--regen-dress": ""}
    """
}

process sampler {
    publishDir "$baseDir/trajs/iter$iter/seed1", mode: 'link'
    label 'pinn'

    input:
    tuple val(iter), path(init), path(model)  from sampleInp

    output:
    tuple val(iter), file('ref.xyz') into trajs

    script:
    """
    #!/usr/bin/env python3
    import pinn
    import numpy as np
    import tensorflow as tf
    from ase import units
    from ase.io import read, write
    from ase.io.trajectory import Trajectory
    from ase.md import MDLogger
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md.nptberendsen import NPTBerendsen

    calc = pinn.get_calc("$model/params.yml")
    for seed in range($params.sampleSeeds):
        rng = np.random.default_rng(seed)
        atoms = read("$init")
        atoms.set_calculator(calc)
        MaxwellBoltzmannDistribution(atoms, 330.*units.kB, rng=rng)
        dt = 0.5 * units.fs
        steps = int($params.sampleTime*1e3*units.fs/dt)
        dyn = NPTBerendsen(atoms, timestep=dt, temperature=330, pressure=1,
                          taut=dt * 100, taup=dt * 1000, compressibility=4.57e-5)
        interval = int($params.sampleInterv*1e3*units.fs/dt)
        dyn.attach(MDLogger(dyn, atoms, 'aug.log', mode="a"), interval=interval)
        dyn.attach(Trajectory('aug.traj', 'a', atoms).write, interval=interval)
        try:
            dyn.run(steps)
        except:
            pass
    traj = read('aug.traj', index=':')
    [atoms.wrap() for atoms in traj]
    write('ref.xyz', traj)
    """
}

process pinnlabel {
    publishDir "$baseDir/trajs/iter$iter/seed$seed", mode: 'link'
    label 'pinn'

    input:
    tuple val(iter), val(seed), path(traj), path(model) from qbcInp

    output:
    tuple val(iter), file('pred.xyz') into qbcOther

    script:
    """
    #!/usr/bin/env python3
    import pinn, yaml
    import tensorflow as tf
    from ase import units
    from ase.io import read, write
    with open('$model/params.yml') as f:
        params = yaml.safe_load(f)
        params['model_dir'] = '$model'
    calc = pinn.get_calc(params)
    traj = read("$traj", index=':')
    with open('pred.xyz', 'w') as f:
        for atoms in traj:
            atoms.wrap()
            atoms.set_calculator(calc)
            atoms.get_potential_energy()
            write(f, atoms, format='extxyz', append='True')
    """
}

process filter {
    publishDir "$baseDir/trajs/iter$iter", mode: 'link'
    label 'pinn'

    input:
    tuple val(iter), path('ds??/*') from qbcRef.mix(qbcOther).groupTuple(size:params.modelSeeds)

    output:
    tuple val{iter}, path('qbc.xyz') into qbcDs
    path('var.xyz')

    script:
    """
    tips qbc ds*/* -o var -of xyz
    tips filter var.xyz $params.qbcFilter -o qbc -of xyz
    """
}

process lmplabel {
    publishDir "$baseDir/datasets/", pattern: "{train,test}*.xyz", mode: 'link'
    publishDir "$baseDir/trajs/iter$iter", pattern: "label.*", mode: 'link'
    label 'lammps'

    input:
    tuple val(iter), path(ds) from qbcDs
    path labeller from Channel.fromPath(params.labeller).collect()

    output:
    tuple val(iter), path("train_${iter}.xyz") into aug4combine
    tuple val{iter}, path('restart.xyz') into restart
    path "test_${iter}.xyz"
    path "label.*"

    script:
    """
    tips convert $ds -o filtered -of 'lammps' --emap '1:1,8:2,11:3,17:4'
    mpirun -np $task.cpus lmp_mpi -in label.lmp || echo LAMMPS aborted
    sed -i '/WARNING/d' label.log
    tips filter label.dump --log label.log -o label -of xyz\
        --emap '1:1,2:8,3:11,4:17' --units real $params.labelFilter
    tips filter label.xyz -o tmp -of xyz $params.restartFilter
    tac tmp.xyz | grep -m1 -A1 -B999 Lattice | tac > restart.xyz
    tips split label.xyz -s 'train_$iter:0.9,test_$iter:0.1' -of xyz
    """
}

process finalmd {
    publishDir "$baseDir/trajs/final", mode: 'link'
    label 'pinn'

    input:
    path model from finalModel.last().map{it[2]}

    output:
    path 'final*'

    script:
    """
    #!/usr/bin/env python3
    import pinn
    import numpy as np
    import tensorflow as tf
    from ase import units
    from ase.io import read, write
    from ase.io.trajectory import Trajectory
    from ase.md import MDLogger
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md.nptberendsen import NPTBerendsen

    calc = pinn.get_calc("$model/params.yml")
    atoms = read("${file(params.sampleInit)}")
    atoms.set_calculator(calc)
    MaxwellBoltzmannDistribution(atoms, 330.*units.kB)
    dt = 0.5 * units.fs
    steps = int(100*1e3*units.fs/dt)
    dyn = NPTBerendsen(atoms, timestep=dt, temperature=330, pressure=1,
                      taut=dt * 100, taup=dt * 1000, compressibility=4.57e-5)
    interval = int(0.1*1e3*units.fs/dt)
    dyn.attach(MDLogger(dyn, atoms, 'final.log', mode="a"), interval=interval)
    dyn.attach(Trajectory('final.traj', 'a', atoms).write, interval=interval)
    try:
        dyn.run(steps)
    except:
        pass
    traj = read('final.traj', index=':')
    [atoms.wrap() for atoms in traj]
    write('final.xyz', traj)
    """
}
