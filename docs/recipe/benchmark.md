# Benchmark workflows

## Usage

The below commands generate a workflow in `main.nf`and then runs it via
nextflow. The TIPS wizard generates a generic skeleton for designing the
workflows, to futher tweak it, please refer to the annotated `main.nf` and
`nextflow.config` files.

=== "commands"

    ```bash
    tips wizard benchmark
    nextflow run main.nf
    ```

=== "main.nf"

    ```groovy
    #!/usr/bin/env nextflow
    
    params.dataset = 'qm9'
    params.input = './inputs/*.yml'
    params.md_init = 'h2o.xyz'
    params.md_flag = '--nvt --T 373 --step 1000 --dt 0.5'
    params.rdf_flag = '--tags OO,OH --rc 5'

    include {pinnTrain} from './pinn.nf'
    include {aseMD} from './ase.nf'
    include {rdf} from './analysis.nf'
    
    workflow {
      dataset = Channel.fromPath(params.dataset)
      input = Channel.fromPath(params.input)
      
      models = pinnTrain(dataset, input)
      trajs = aseMD(models)
      rdf(trajs, rdf_tags)
    }
    ```

=== "nextflow.config"

    ```groovy
    profiles {
      standard {
        process {
          cpus=1
          errorStrategy='ignore'
          withLabel: pinn {container='yqshao/pinn:master'}
          withLabel: tame {container='yqshao/tame:master'}
        }
        executor {
          name = 'local'
          cpus = 16
        }
      }
    ```

## Other links

- [List of available datasets]()
- [Sample model input files]()
- [Adding a custom analysis]()
