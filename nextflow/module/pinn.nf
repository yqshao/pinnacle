nextflow.enable.dsl=2

params.publish = 'pinn'

process train {
  label 'pinn'
  publishDir "$params.publish/$name"

  input:
    tuple val(name), path(dataset), path(input, stageAs:'input'), val(flags)

  output:
    tuple val(name), path('model', type:'dir'), emit: model
    tuple val(name), path('pinn.log'), emit: log

  script:
    convert_flag = "${(flags =~ /--seed[\s,\=]\d+/)[0]}"
    train_flags = "${flags.replaceAll(/--seed[\s,\=]\d+/, '')}"
    dataset = (dataset instanceof Path) ? dataset : dataset[0].baseName+'.yml'
    """
    #!/bin/bash

    pinn convert $dataset -o 'train:9,eval:1' $convert_flag

    if [ ! -f $input/params.yml ];  then
        mkdir -p model; cp $input model/params.yml
    else
        cp -rL $input model
    fi
    pinn train model/params.yml --model-dir='model'\
        --train-ds='train.yml' --eval-ds='eval.yml'\
        $train_flags
    pinn log model/eval > pinn.log
    """
}

process md {
  label 'pinn'
  publishDir "$params.publish/$name"

  input:
    tuple val(name), path(model,stageAs:'model*'), path(init, stageAs:'init*'), val(flags)

  output:
    tuple val(name), path('asemd.traj'), emit: traj
    tuple val(name), path('asemd.log'), emit: log

  script:
    """
    #!/usr/bin/env python
    import re
    import pinn
    import tensorflow as tf
    from ase import units
    from ase.io import read
    from ase.io.trajectory import Trajectory
    from ase.md import MDLogger
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md.nptberendsen import NPTBerendsen
    from ase.md.nvtberendsen import NVTBerendsen
    from tips.bias import EnsembleBiasedCalculator

    # ------------ patch ase properties to write extra cols --------------------
    from ase.calculators.calculator import all_properties
    all_properties+=[f'{prop}_{extra}' for prop in ['energy', 'forces', 'stress'] for extra in ['avg','std','bias']]
    # --------------------------------------------------------------------------

    setup = {
      'ensemble': 'nvt', # ensemble
      'T': 340, # temperature in K
      't': 50, # time in ps
      'dt': 0.5, # timestep is fs
      'taut': 100, # thermostat damping in steps
      'taup': 1000, # barastat dampling in steps
      'log-every': 20, # log interval in steps
      'pressure': 1, # pressure in bar
      'compressibility': 4.57e-4, # compressibility in bar^{-1}
      'bias': None,
      'kb': 0,
      'sigma0': 0,
    }

    flags = {
      k: v for k,v in
        re.findall('--(.*?)[\\s,\\=]([^\\s]*)', "$flags")
    }
    setup.update(flags)
    ensemble=setup['ensemble']
    T=float(setup['T'])
    t=float(setup['t'])*units.fs*1e3
    dt=float(setup['dt'])*units.fs
    taut=int(setup['taut'])
    taup=int(setup['taup'])
    every=int(setup['log-every'])
    pressure=float(setup['pressure'])
    compressibility=float(setup['compressibility'])

    ${(model instanceof Path) ?
    "calc = pinn.get_calc('$model')" :
    """
    models = ["${model.join('", "')}"]
    calcs = [pinn.get_calc(model) for model in models]
    if len(calcs) == 1:
        calc =  calcs[0]
    else:
        calc = EnsembleBiasedCalculator(calcs,
                                        bias=setup['bias'],
                                        kb=float(setup['kb']),
                                        sigma0=float(setup['sigma0']))
    """}

    atoms = read("$init")
    atoms.set_calculator(calc)
    if not atoms.has('momenta'):
        MaxwellBoltzmannDistribution(atoms, T*units.kB)

    if ensemble == 'npt':
        dyn = NPTBerendsen(atoms, timestep=dt, temperature=T, pressure=pressure,
                      taut=dt * taut, taup=dt * taup, compressibility=compressibility)
    if ensemble == 'nvt':
        dyn = NVTBerendsen(atoms, timestep=dt, temperature=T, taut=dt * taut)

    dyn.attach(
        MDLogger(dyn, atoms, 'asemd.log',stress=True, mode="w"),
        interval=int(every))
    dyn.attach(
        Trajectory('asemd.traj', 'w', atoms).write,
        interval=int(every))
    dyn.run(int(t/dt))
    """
}
