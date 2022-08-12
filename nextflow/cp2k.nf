nextflow.enable.dsl=2

params.publish = 'cp2k'
params.cp2k_cmd = 'cp2k'
params.cp2k_aux = null

process cp2k {
  tag "$name"
  label 'cp2k'
  publishDir "$params.publish/$name"

  input:
    tuple val(name), path(input), path(aux)

  output:
    tuple val(name), path('*.{ener,xyz,stress}'), emit:traj
    tuple val(name), path('cp2k.log'), emit:logs
    tuple val(name), path('*.restart'), emit:restart, optional:true

  script:
    """
    #!/bin/bash
    $params.cp2k_cmd -i $input | tee cp2k.log
    """
}


process cp2kGenInp {
  tag "$name"
  label 'tips'
  publishDir "$params.publish/$name"

  input:
    tuple val(name), path(input,stageAs: 'cp2k_skel.inp'), path(init), val(flags)

  output:
    tuple val(name), path('cp2k.inp')

  script:
  """
  #!/usr/bin/env python
  import re
  from ase.data import chemical_symbols as symbol
  from tips.io import load_ds
  # read flags
  setup = {
    'emap': None,
    'idx': -1,
    'fmt': 'auto',
  }
  flags = {
    k: v for k,v in
      re.findall('--(.*?)[\\s,\\=]([^\\s]*)', "$flags")
  }
  setup.update(flags)
  # load the desired geometry
  ds = load_ds("$init", fmt=flags['fmt'])
  if setup['emap'] is not None:
      ds = ds.map_elems(setup['emap'])
  datum = ds[int(setup['idx'])]
  # edit the input file
  coord = [f'  {symbol[e]} {x} {y} {z}' for e, (x,y,z) in zip(datum['elem'], datum['coord'])]
  cell = [f'  {v} {x} {y} {z}' for v, (x, y, z) in zip('ABC', datum['cell'])]
  subsys = ['&COORD']+coord+['&END COORD']+['&CELL']+cell+['&END CELL']
  lines = open("$input").readlines()
  for idx, line in enumerate(lines):
      if '&END SUBSYS' in line:
          indent = len(line) - len(line.lstrip())
          break
  subsys = [' '*(indent+2) + l + '\\n' for l in subsys]
  lines = lines[:idx] + subsys + lines[idx:]
  with open('cp2k.inp', 'w') as f:
      f.writelines(lines)
  """
}


workflow cp2kMD {
  take:
    ch // [name, input, init, flags]

  main:
    ch | cp2kGenInp // -> [name, inp]
       | map {name, inp -> [name, inp, file(params.cp2k_aux)]}
       | cp2k

  emit:
    traj = cp2k.out.traj
    logs = cp2k.out.logs
    restart = cp2k.out.restart
}
