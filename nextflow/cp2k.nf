nextflow.enable.dsl=2

params.cp2k_cmd = 'cp2k'
params.cp2k_aux = null

process cp2k {
  label 'cp2k'
  publishDir "$publish"

  input:
    val name
    path input
    path init
    path aux
    val publish

  output:
    tuple val(name), path('*.{ener,xyz,stress}'), emit:traj
    tuple val(name), path('cp2k.log'), emit:log
    tuple val(name), path('*.restart'), emit:restart, optional:true

  script:
  """
  #!/bin/bash
  $params.cp2k_cmd -i $input | tee cp2k.log
  """
}


process ck2pInp {
  label 'tips'
  cache false
  publishDir "$publish"

  input:
    val name
    patth input
    path init

  output:
    tuple val(name), path('generated.inp')

  script:
  """
  #!/bin/bash

  """
}


workflow cp2kMD {
  take:
    name
    input
    init
    publish

  main:
    ch_inp = cp2kInp(name, input, init)
    ch_pub = name.merget(publish)
    ch_inp
      .join(ch_pub)
      .multiMap{
        name, inp, publish ->
        name: name
        inp: inp
        aux: file(params.cp2k_aux)
        publish: publish
      }
      .set {ch}
    out = cp2k(ch.name, ch.inp, ch.aux, ch.publish)

  emit:
    traj = out.traj
    log = out.log
    restart = out.restart
}
