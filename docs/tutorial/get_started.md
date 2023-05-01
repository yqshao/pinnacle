# Get Started

In this demo, we prepare a dataset with the xTB for different water/ethanol
mixtures, and run an activated learning workflow starting from the initial
dataset.

## Installation

PiNNAcLe depends on nextflow, it can be installed following the official
[documentation][nf_install]. By default, the workflows will be run through
singularity containers, and you need to install that as well, see their
[documentation][sig_install].

Singularity is the recommended way of running the workflows, it is intended for
running and sharing software, and avoids complications with compilation and
environment setup. But if you prefer otherwise, you can easily configure
workflow to run with a different profile, which controls how and on what
resource each process is run, read more in the [Profiles] section.

[nf_install]: https://www.nextflow.io/docs/latest/getstarted.html#installation
[sig_install]: https://docs.sylabs.io/guides/main/user-guide/quick_start.html#quick-installation-steps
[Profiles]: ../../profiles/overview

## Get the workflow

The PiNNAcLe workflow can be downloaded as a [copier] template. The template
simply downloads the PiNNAcLe modules and sets up a project directory ready to
be extended. If you do not have copier installed, you can simply clone the repo.

[copier]: https://copier.readthedocs.io/en/stable/#installation

```bash
copier gh:teoroo-cmc/pinnacle demo_proj
# or git clone https://github.com/Teoroo-CMC/PiNNAcLe.git demo_proj
```

You may also run the workflow directly as a shared pipeline, as demonstrated
[here][pipeline], we chose to have the modules explicitly in the project folder
so that you can see how they can be adapted to your need in later sections of
the tutorial.

[pipeline]: https://www.nextflow.io/docs/latest/sharing.html#running-a-pipeline

## Running the workflow

To run the workflow:

```bash
nextflow run main.nf # add -profile your_profile if necessary
```

You might need to wait a while for the singularity images to be pulled,
afterward you should see something like this:

```
Launching `main.nf` [friendly_carlsson] DSL2 - revision: e28d0f3d5d
executor >  local (18)
[0d/2b1c2a] process > we_demo:mol2box (2)        [100%] 9 of 9 ✔
[3b/ebc632] process > we_demo:md:aseMD (3)       [  0%] 0 of 6
[-        ] process > we_demo:mkds               -
[-        ] process > we_demo:mkgeo              -
[-        ] process > we_demo:acle:loop:train    -
[-        ] process > we_demo:acle:loop:md       -
[-        ] process > we_demo:acle:loop:convert  -
[-        ] process > we_demo:acle:loop:sp:aseSP -
[-        ] process > we_demo:acle:loop:merge    -
[-        ] process > we_demo:acle:loop:check    -
[-        ] process > we_demo:acle:loop:dsmix    -
```

What listed above as "processes" are the basic element of nextflow workflows. A
workflow defines how different processes are linked, and the configuration
controls how they are launched in practice (e.g., environment setup,
parallelization, resource allocation, etc.)

You might notice that the processes are tracked by a 2-level hex code such as
`0d/2d1c2a`. This is how nextflow arrange the output of processes, each process
will be launched in a folder named after this unique code. The code will be
useful for restarting the workflow and for debugging the processes, the folder
structure is illustrated below.


```bash
demo_proj
├── main.nf # main entry point for the workflows
├── nextflow # nextflow recipes and modules
│   ├── acle.nf
│   ├── we_demo.nf
│   └── module
│       ├── pinn.nf
│       ├── tips.nf
│       └── ...
├── we_demo # output folders are linke to work directory
│   ├── init
│   │   ├── geo
│   │   │   ├── e32-1.xyz -> /xx/demo_proj/work/c6/322d1c1...
│   │   └── ...
│   └── ...
└── work # actual execution directory of the workflow
    ├── c6
    │   ├── 322d1c1f55ec381f34c361bc4e743f
    │   │   ├── e32-1.xyz
    │   │   ├── .commnad.sh
    │   │   ├── .commnad.run
    │   │   ├── .commnad.out
    │   │   └── ...
    │   └── ...
    └── ...
```


## Useful commands

Before going further, below are some nextflow commands that you might find
useful. For a more comprehensive description, check the nextflow [CLI
reference].

[CLI reference]: https://www.nextflow.io/docs/latest/cli.html


=== "Checking"

    Nextflow keeps track of what you have run and how your run them
    Each run of the workflow is attached to a "run name", with which you can resume, log,
    and clean up the jobs, try `nextflow log tiny_brattain` (replace with your own run) to
    get a list of all jobs in a run.

    ```
    $ nextflow log 
    ---
    TIMESTAMP              DURATION        RUN NAME                    STATUS  REVISION ID     SESSION ID                              COMMAND
    2022-04-26 22:52:18    1h 40m 53s      tiny_brattain               OK      e7132a82d7      bf527311-a5d9-4f7b-b615-d50ff99e6ec5    nextflow run main.nf
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
    [85/e10f9e] process > train (1) [50%] 3 of 6, cached: 3
    ```

=== "Cleaning"

    Keeping a complete record of things has a cost, the script and 
    intermediate results for past runs might accumulate to a large number of file. 
    It is recommended to clean up your unneeded runs regularly. If you are 
    restricted by disk, you should also consider setting the [scratch]
    directive.
    
    [scratch]: https://www.nextflow.io/docs/latest/process.html#scratch

    ```
    $ nextflow clean -f disturbed_crick
    ---
    Removed /home/user/proj/work/3f/f350726963c2372d1673b900553d6f
    Removed /home/user/proj/work/83/60f47f79880256257fc53a19b6ffac
    Removed /home/user/proj/work/5f/bd7e2a391cb53a4e20057c4641e408
    ```

## What next?

- [Workflow](../workflow): more on how to run workflows;
- [Configure](../configure): how the processes are executed.
