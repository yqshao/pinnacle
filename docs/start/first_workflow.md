# Your first workflow

In this section, you'll run the workflow you obtain from the previous step and
learn some basics about

## Running the workflow

To the workflow you get from the preview section, simply run:

=== "Command"

    ```
    nextflow run main.nf
    ```

if you have chosen the default config, nextflow will fetch the necessary
singularity images and run the training on you local computer. If you are lucky,
you will see the trained models as well as their training logs in you `models`
folder.

## Structure of main.nf

To see what is happening here, we first look at the `main.nf` file:

=== "main.nf"

    ```groovy
    #!/bin/env nextflow

    // Nextflow script for training and evaluation MLPs,
    // generated with TIPS v0.1.0 at 21:32-220503

    import {trainer} from 'tips/pinn.nf'
    import {convert} from 'tips/tips.nf'

    params.ds = 'datasets/qm9' // input data
    params.inps = 'inputs/pinet.yml' // model inputs
    params.seeds = '1,2,3,4,5' // random sees to initialze
    params.splits = '90:10' // train,eval splits

    workflow {
      ds = Channel.fromPath(params.ds)
      input = Channel.fromPath(params.inps)
      splits = Channel.fromList(params.split)

      names = ds.combine(inputs).combine(splits).combine(seeds)
        .map{ds, inp, split, seed -> "$ds-$inp.name-$split-$seed"}
      dataset = convert(ds, params.splits)
      trainer(dataset, input, names)
    }
    ```

This file defines the workflow, as well as the adjustable parameters. For
instance, the "ds" parameter can be adjusted at runtime with `nextflow run --ds datasets/water.yml`.

You might notice that this workflow imports some "processes" (such as trainer)
from the other "modules". Looking into `pinn.nf`, you'll find that it defines
the actual commands ran to perform the training, and the input/output thereof.

In TIPS, the processes can be used interchangably so long as they share a
similar input/output pattern. For example, you can easily integrate a
Python-based training script into any existing workflow, so long as it consumes
a dataset and output a model.

## Nextflow basics

To list previous runs in the project folder:

=== "Command"

    ```
    nextflow log
    ```

=== "Output"

    ```
     TIMESTAMP              DURATION        RUN NAME                    STATUS  REVISION ID     SESSION ID                              COMMAND
     2022-04-26 22:28:35    1m 13s          tender_varahamihira         OK      e7132a82d7      4f445e64-8d54-48dd-9fea-b57d9be3e5c9    nextflow run
     2022-04-26 22:52:18    1h 40m 53s      tiny_brattain               OK      e7132a82d7      bf527311-a5d9-4f7b-b615-d50ff99e6ec5    nextflow run
    ```

To restart from a past run with new parameters:

=== "Command"

    ```
    nextflow run --inps 'inputs/*.yml' -resume
    ```

=== "Output"

    ```
    Launching `main.nf` [disturbed_crick] - revision: 3439ecc683
    executor > local (6)
    [85/e10f9e] process > trainner (3) [50%] 3 of 6, cached: 3
    ```

Note that nextflow automatically caches your completed tasks, with the `-resume`
command, the tasks with the exact same input will be resued. For a more
comprehensive description, check the Nextflow CLI
[reference](https://www.nextflow.io/docs/latest/cli.html).

## What next?

You might notice that we did not touch upon any thing related to the
environment, e.g., what if I need a different version of PiNN, or if I with to
run the jobs in a queuing system? This is intentional since in TIPS the workflow
is decoupled from the "configuration" of computational environment as much as
possible. This design choice is to minimize the hassle when using the workflow.
In the next section, you will get familiar with the `nextflow.config` file,
where the computational resource, environment will be defined.
