nextflow.enable.dsl=2

params.publish = ''
params.ase_aux = '/dev/null'
params.ase_label = 'ase'
params.ase_props = '["energy", "forces"]'

workflow md {
  take:
    ch // [name, inp, geo, flags]

  main:
    ch | map {name, inp, geo, flags -> \
              [name, inp, geo, flags, file(params.ase_aux)]} \
       | aseMD

  emit:
    traj = aseMD.out.traj
    logs = aseMD.out.logs
}

workflow sp {
  take:
    ch // [name, geo]

  main:
    ch | map {name, inp, geo -> \
              [name, inp, geo, file(params.ase_aux) ]} \
       | aseSP

  emit:
    aseSP.out
}


// --8<-- [start:sp]
process aseSP {
  label "$params.ase_label"
  publishDir "$params.publish/$name"

  input:
    tuple val(name), path(calc, stageAs: 'calc.py'), path(geo), path(aux)

  output:
    tuple val(name), path("sp.xyz")

  script:
    """
    #!/usr/bin/env python
    from ase.io import read, write
    from calc import calc

    atoms = read('$geo')
    atoms.calc=calc
    calc.calculate(atoms, properties=${params.ase_props})

    write('sp.xyz', atoms)
    """
}
// --8<-- [end:sp]


process aseMD {
  label "$params.ase_label"
  publishDir "$params.publish/$name"

  input:
    tuple val(name), path(calc, stageAs:'calc.py'), path(init), val(flags), path(aux)

  output:
    tuple val(name), path('asemd.traj'), emit: traj
    tuple val(name), path('asemd.log'), emit: logs

  script:
    """
    #!/usr/bin/env python
    from calc import calc
    import re
    from ase import units
    from ase.io import read
    from ase.io.trajectory import Trajectory
    from ase.md import MDLogger
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md.nptberendsen import NPTBerendsen
    from ase.md.nvtberendsen import NVTBerendsen

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
