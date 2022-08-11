nextflow.enable.dsl=2

params.lmp_cmd = 'lmp'

process lammpsMD {
  tag "$name"
  publishDir "$publish"
  label 'lammps'

  input:
    val name
    path input
    path aux
    val publish

  output:
    tuple val(name), path('*.dump'), emit: traj
    tuple val(name), path('log.lammps'), emit: log

  script:
    """
    #!/bin/bash
    $params.lmp_cmd -i $input
    """
}
