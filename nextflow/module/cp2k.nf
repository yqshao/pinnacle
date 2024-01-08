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
    tuple val(name), path('cp2k.log'), emit:logs
    tuple val(name), path('*.{ener,xyz,stress,cell}'), emit:traj, optional: true
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
    tuple val(name), path(input,stageAs: 'cp2k_skel.inp'), path(ds), val(flags)

  output:
    tuple val(name), path('*.inp')

  script:
  """
  tips cp2kinp $input $ds $flags
  """
}


workflow md {
  take:
    ch // [name, input, init, flags]

  main:
    ch | cp2kGenInp // -> [name, inp] \
       | map {name, inp -> [name, inp, file(params.cp2k_aux)]} \
       | cp2k

  emit:
    traj = cp2k.out.traj
    logs = cp2k.out.logs
    restart = cp2k.out.restart
}

workflow sp {
  take:
    ch // [name, inp, geo]

  main:
    ch | map {name, inp, geo -> \
              [name, inp, geo, '-f asetraj']} \
       | cp2kGenInp \
       | map {name, inp -> [name, inp, file(params.cp2k_aux)]} \
       | cp2k

  emit:
    cp2k.out.logs
}
