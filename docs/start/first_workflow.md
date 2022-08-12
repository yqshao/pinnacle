# Your first workflow

In this section, you'll run the workflow from the previous step and see what
TIPS can do.

## Running the workflow

To the workflow you get from the preview section, simply run:

```bash
nextflow run main.nf
```

if you have chosen the default config, nextflow will fetch the necessary
singularity images and run the training on you local computer. If you are lucky,
you will see the trained models as well as their training logs in you `models`
folder.

## Structure of the project

The TIPS script provides a way to run your script in a reproducible, and
scalable way. With the `main.nf` script, you are ready to run the same workflow
on any computational resource, and easily include your own tricks to fit your
need. See the annotated fodler structure and files below for more explanation.

=== "Folder structure"

    ```bash
    your_
    |- main.nf # (1)
    |- nextflow.config # (2)
    |- tips # (3)
       |- pinn.nf 
    |- work # (4)
       |- f3/fa4ce305c2b0d057af8b83ad9b0aec
       |- ....
    |- models # (5)
       |- qm9-pinet-1
       |- qm9-pinet-2
       |- ....
    ```

    1. This is the workflow file that you execute, see the tab on the left.
    2. This controls how the jobs are run on different resources.
    3. TIPS provides several `modules` for different tasks, those are stored in the `tips` subfolder,
       see the [modules documentation](../recipe/overview) for their specific usage. 
    4. This holds the actual runtime folder of the jobs, where the scripts, data and intermediate results sites.
    5. The `models` (or alike) folders are "published" by nextflow, tips have a sysmtatic way to order the jobs, see [general rules](workflow/#general_rules).

=== "main.nf"

    ```groovy
    import { pinnTrain } from 'tips/pinn.nf' // (1)


    params.ds = 'datasets/h2o.data' // (2)
    params.inps = 'inputs/pinet.yml'
    params.flags = '--train-steps 1000000 --batch 30' 

    workflow {
      channel.fromPath(params.ds) 
        .combine(Channel.fromPath(params.inps))
        .multiMap {
          ds, inp:
          name: $"$ds.baseName-inp.baseName"
          ds: ds
          inp: inp
          flags: params.flags
        }
        .set { ch } // (3)
      pinnTrain (ch.name, ch.ds, ch.inp, ch.flags)
    }
    ```

    1. This imports the `pinnTrain` process that does the actual training, you can
        find from the [documentation](../module/pinn/#pinntrain) that under the 
        hood, it's just a bash script. To build your own workflow, feel free to
        substitue in your favorate ML/MD package.
    2. The `params.xx` lines defines variables you can change during runtime, for
       instance running the same script with `nextflow run main --ds 'dataset/water.data'` 
       trains your model on a different dataset.
    3. TIPS build up complex workflows following a few rules. Here, we create the data channels
       and name the runs according to the dataset and PiNN input. The channels are built from named channels with the multiMap function, see [nextflow
       documentation](https://www.nextflow.io/docs/latest/operator.html#multimap)
       for the details.

## Nextflow basics

Below going further, there are some nextflow tricks that you might find useful:

=== "Checking"

    ``` bash
    $ nextflow log # (1)
    ------------------------------------ # (2)
    TIMESTAMP              DURATION        RUN NAME                    STATUS  REVISION ID     SESSION ID                              COMMAND
    2022-04-26 22:28:35    1m 13s          tender_varahamihira         OK      e7132a82d7      4f445e64-8d54-48dd-9fea-b57d9be3e5c9    nextflow run
    2022-04-26 22:52:18    1h 40m 53s      tiny_brattain               OK      e7132a82d7      bf527311-a5d9-4f7b-b615-d50ff99e6ec5    nextflow run
    ```

    1. Nextflow keeps track of what you have run and how your run them; to check the record, use `nextflow log`.
    2. Each run of the workflow is attached to a "run name", with which you can resume, log, 
       and clean up the jobs, try `nextflow log tiny_brattain` (replace with your own run) to 
       get a list of all jobs in a run.

=== "Resuming"

    ```bash
    $ nextflow run main.nf -resume # (1)
    ---
    Launching `main.nf` [disturbed_crick] - revision: 3439ecc683
    executor > local (6)
    [85/e10f9e] process > pinnTrain (h2o-pinet-1) [50%] 3 of 6, cached: 3
    ```
    
    1. Nextflow caches the jobs according to their inputs, if your run is interrupted, you can resume from where your stopped with the `-resume`. You can also specify a previous run 
       to resume from: `nextflow run main.nf -resume tiny_brattain`

=== "Cleaning"

    ```bash
    $ nextflow clean -f disturbed_crick # (1)
    --- 
    Removed /home/user/proj/work/3f/f350726963c2372d1673b900553d6f
    Removed /home/user/proj/work/83/60f47f79880256257fc53a19b6ffac
    Removed /home/user/proj/work/5f/bd7e2a391cb53a4e20057c4641e408
    ```
    
    1. Keeping a complete record of things come at a cost, you are 
       recommended to clean up your unneeded runs regularly. There are also 
       the `-before`, `-after`and `-but` to clean up multiple runs a time.

For a more comprehensive description, check the [Nextflow CLI
reference](https://www.nextflow.io/docs/latest/cli.html).

## What next?

- [Configure your run](configure) your run to work on any computational
  resource.
- Some more complex examples of how [workflows](../../recipe/overview) are
  composed.
- Explore how TIPS as a python [library](../../library).
