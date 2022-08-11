nextflow.enable.dsl=2

process convert {
  label 'tips'
  publishDir "$publish"

  input:
    val name
    path input
    val flags
    val publish

  output:
    tuple val(name), path('converted.*')

  script:
    """
    #!/bin/bash
    tips covnert $input $flags
    """
}
