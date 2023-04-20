nextflow.enable.dsl=2

params.publish = "convert"

process convert {
  label 'tips'
  publishDir "$params.publish/$name"

  input:
    tuple val(name), path(input), val(flags)

  output:
    tuple val(name), path('converted.*')

  script:
    """
    #!/bin/bash
    tips covnert $input $flags
    """
}
