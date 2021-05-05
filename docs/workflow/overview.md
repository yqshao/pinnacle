# Workflow Overview

TIPS provides a number of workflows to perform atomistic machine learning (AML)
tasks. Those include:

- Strategies (e.g. `explore`) contains specific AML workflows like active
  learning, delta ML, etc.
- Implementations (e.g. `pinn`, `lammps`) contains implementation of workflows.
- The `adaptor` provides abstract workflows that are implementation agnostic.

## Usage

The aforementioned workflows are written in the Nextflow language, with which
the workflows can be automatically pulled from Github using the following
command:

```shell
nf run yqshao/tips --mode explore --initDs myDs.data
```

The `--mode` options specifies what workflow to run, it can either be a strategy
or an abstract workflow.

## Installation

To modify an existing workflow to import it as part for your own workflow, you
will need to install the TIPS workflows. One way to do that is simply
downloading the `.nf` files to you local computer.

```shell
mkdir -p tips && curl https://codeload.github.com/yqshao/tips/tar.gz/master | \
  tar -C tips -xz --strip=2 tips-master/nextflow/
```

The local workflows can be executed as such:

```shell
nf run tips/explore.nf --initDs myDs.data
```

## Inclusion

One potential purpose of installing TIPS workflows is to use them as a component
of your own workflow, for instance, the following workflow imports the explore
strategy and tries it with different training parameters.

```groovy
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {explore} from './tips/explore'

meta = Channel.value(null)
inputs = Channel.of('pinet.yml', 'bpnn.yml')
                .map{[trainInp: it, subDir: "$it.baseName"]}

outputs = explore(meta.combine(inputs))
```

All workflows and processes in TIPS follows a same pattern for input/output, all
inputs should be a tuple of `[meta, inputs]`, and outputs `[meta, output]`. The
metadata are simply copied to the output, mainly for the outputs to be
identified by the outer workflow. 

Most options of the a workflow, like the `trainInp` here, can be adjusted either
through the CLI `params`, or through the its input channel. You might want to
check the documentation for the specific workflow for the available otpions.

