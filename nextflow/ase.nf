nextflow.enable.dsl=2

params.publish = 'ase'

process aseMD {
  tag "$name"
  label 'ase'
  publishDir "$param.publish/$name"

  input:
    tuple val(name), path(model,stageAs:'model*'), path(init), val(flags)

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
    from ase.calculators.mixing import AverageCalculator
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md.nptberendsen import NPTBerendsen

    setup = {
      'ensemble': 'npt', # ensemble
      'T': 373, # temperature in K
      't': 100, # time in ps
      'dt': 0.5, # timestep is fs
      'taut': 100, # thermostat damping in steps
      'taup': 1000, # barastat dampling in steps
      'log-every': 5, # log interval in steps
      'pressure': 1, # pressure in bar
      'compressibility': 4.57e-4 # compressibility in bar^{-1}
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
    calc = AverageCalculator(calcs)
    """}

    atoms = read("$init")
    atoms.set_calculator(calc)
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
