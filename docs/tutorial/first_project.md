# Your first PiNNAcLe workflow

In this section, you will create a project with PiNNAcLe and see how it works.

## Using a template

```bash
copier gh:teoroo-cmc/pinnacle
```

## Running the workflow

To the workflow you get from the preview section, simply run:

```bash
nextflow run main.nf
```

If you have chosen the default config, nextflow will fetch the necessary
singularity images and run the training on you local computer.

## Structure of the project

The PiNNAcLe templates arranges your workflow in a reproducible, and scalable
way. With the `main.nf` script, you are ready to run the same workflow on any
computational resource, and easily include new modules. See the annotated folder
structure and files below for more explanation.

=== "Folder structure"

    ```bash
    my_proj
    |- main.nf # This is the workflow file that you execute, see the tab on the left.
    |- nextflow.config # This controls how the jobs are run on different resources.
    |- nf-module # this contains the nextflow modules that PiNNAcLe provides. 
       |- pinn.nf
    |- work # This holds the work folder of the jobs, along with intermediate results
       |- f3/fa4ce305c2b0d057af8b83ad9b0aec
       |- ....
    |- models # Those are the outputs "published" by nextflow, linked to their work folders
       |- qm9-pinet-1
       |- qm9-pinet-2
       |- ....
    ```

=== "main.nf"

    ```groovy
    import { pinnTrain } from 'nf-module/pinn.nf'

    // The `params.xx` lines defines variables you can change during runtime, for
    // instance running the same script with `nextflow run main --ds 'dataset/water.data'`
    // trains your model on a different dataset.
    params.ds = 'datasets/h2o.data'
    params.inps = 'inputs/pinet.yml'
    params.flags = '--train-steps 1000000 --batch 30'

    workflow {
      channel.fromPath(params.ds) \
        | combine(channel.fromPath(params.inps)) \
        | map {ds, inp -> ["$ds.baseName-$inp.baseName", ds, inp, params.flags]} \
        | pinnTrain
    }
    // PiNNAcLe build up complex workflows following a few rules. You see above that we use
    // a map operator to convert the inputs into a `tuple` of the inputs for `pinnTrain`. 
    // In most cases, we use first item in the tuple as the identifier of jobs.
    ```


## Nextflow basics

Before going further, below are some nextflow commands that you might find
useful. For a more comprehensive description, check the [nextflow CLI
reference](https://www.nextflow.io/docs/latest/cli.html).


=== "Checking"

    Nextflow keeps track of what you have run and how your run them
    Each run of the workflow is attached to a "run name", with which you can resume, log,
    and clean up the jobs, try `nextflow log tiny_brattain` (replace with your own run) to
    get a list of all jobs in a run.
    ```
    $ nextflow log 
    ---
    TIMESTAMP              DURATION        RUN NAME                    STATUS  REVISION ID     SESSION ID                              COMMAND
    2022-04-26 22:28:35    1m 13s          tender_varahamihira         OK      e7132a82d7      4f445e64-8d54-48dd-9fea-b57d9be3e5c9    nextflow run
    2022-04-26 22:52:18    1h 40m 53s      tiny_brattain               OK      e7132a82d7      bf527311-a5d9-4f7b-b615-d50ff99e6ec5    nextflow run
    ```

=== "Resuming"

    Nextflow caches the jobs according to their inputs, if your run is
    interrupted, you can resume from where your stopped with the `-resume`. You
    can also specify a previous sessuib to resume from:
    `nextflow run main.nf -resume SESSION_ID`

    ```
    $ nextflow run main.nf -resume
    ---
    Launching `main.nf` [disturbed_crick] - revision: 3439ecc683
    executor > local (6)
    [85/e10f9e] process > pinnTrain (h2o-pinet-1) [50%] 3 of 6, cached: 3
    ```

=== "Cleaning"

    Keeping a complete record of things come at a cost, you are
    recommended to clean up your unneeded runs regularly. There are also
    the `-before`, `-after`and `-but` to clean up multiple runs a time.

    ```
    $ nextflow clean -f disturbed_crick
    ---
    Removed /home/user/proj/work/3f/f350726963c2372d1673b900553d6f
    Removed /home/user/proj/work/83/60f47f79880256257fc53a19b6ffac
    Removed /home/user/proj/work/5f/bd7e2a391cb53a4e20057c4641e408
    ```

## What next?

- [Configure](../configure) your run to work on any computational
  resource;
