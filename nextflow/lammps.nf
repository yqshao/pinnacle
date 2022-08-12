nextflow.enable.dsl=2

params.lmp_cmd = 'lmp'
params.publish = 'lmp'

process lammpsMD {
  tag "$name"
  publishDir "$params.publish/$name"
  label 'lammps'

  input:
    tuple val(name), path(input), path(aux)

  output:
    tuple val(name), path('*.dump'), emit: traj
    tuple val(name), path('log.lammps'), emit: logs

  script:
    """
    #!/bin/bash
    $params.lmp_cmd -i $input
    """
}
