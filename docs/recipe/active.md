# Active workflows

## Usage

The below commands generate a workflow in `main.nf`and then runs it via
nextflow. The TIPS wizard generates a generic skeleton for designing the
workflows, to futher tweak it, please refer to the annotated `main.nf` and
`nextflow.config` files.

=== "commands"

    ```bash
    tips wizard active
    nextflow run main.nf
    ```

=== "main.nf"

    ```groovy
    #!/usr/bin/env nextflow
    
    params.dataset = 'qm9'
    params.input = './inputs/*.yml'
    params.md_init = 'h2o.xyz'
    params.md_flag = '--nvt --T 373 --step 1000 --dt 0.5'

    // to be written
    ```

=== "nextflow.config"

    ```groovygg
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
    
## Strategies

TIPS provides different strategies for an active leraning task, which affects
the efficiency and the resulting NN.

** ENN **: In an ENN (ensemble NN) workflow, the trajectory is propogated with
an ensemble of NN models, while a subset of the trajectory is labelled according
to a given uncertainty tolerance.

** DAS **: The DAS (density-based adaptive sampling) scheme samples the
configuration space by actively biasing the dynamics according to the data
distribution in latent space.

## Other links

- [Sample model input files]()
- [Adding a custom analysis]()
