# Abstract workflows

Three types of abstract workflows are implemented in TIPS, they are workflows
that follows certain input/output patterns, additional parameters can be present
in the inputs if necessary.

| Type     | Required Inputs                                | Output  |
|----------|------------------------------------------------|---------|
| trainer  | `[inp: model, ds:dataset, maxIter: iter, ...]` | models  |
| sampler  | `[inp: model, init: structure, ...]`           | dataset |
| labeller | `[inp: model, ds:dataset, ... ]`               | dataset |
| filter   | `[inp: filter, ds:dataset, ...]`               | dataset |

The structure and dataset files (`ds` or `init`) are converted by the TIPS CLI
when they are parsed across workflows, in principle any file that can be
consumed by the TIPS CLI should work. On the other hand, the `inp` are inputs
that are usually specific to the exact process used.

## Usage

The abstract workflows can be retrieve from the `adaptor` module, where the
exact version retrived according to `param`. For instance, the below script
trains two models with the same datasets with two different programs.

```groovy
params.trainer = 'pinn'
include {trainer as pinn} from './tips/adapter'
include {trainer as runner} from './tips/adapter', addParams(trainer:'runner')

meta = Channel.value(null)
inputs = Channel.of([ds: 'train.xyz'])

workflow{
  inputs | map{it+[subDir: 'pinn']}   | meta.combine | pinn
  inputs | map{it+[subDir: 'runner']} | meta.combine | runner
}
```

## Implemented

Below is a table of the abstract workflows available in TIPS.

| Name   | Processes                              | Description                 |
|--------|----------------------------------------|-----------------------------|
| tips   | `tipsFilter`                           | Filters implemented in TIPS |
| pinn   | `pinnTrain`, `pinnSample`, `pinnLabel` | The PiNN AML package        |
| lammps | `lammpsSample`, `lammpsLabel`          | The LAMMPS MD code          |

Their description and extra options can be found in the implementations section.

## Customization
It is also possible to specify a custom workflow in the workflow by changing the
`trainer`, `sampler` or `labeller` parameter to a relative path starting with
`./`, in this case, adaptor will try to get such a general workflow from the
script, in that case, the workflow or process must be named as `trainer`,
`sampler` or `labeller` accordingly.

```groovy
params.trainer = './custom'
include {trainer} from './tips/adapter'
```

This is useful if you would like to reuse an active learning workflow, but
replace certain module. More information can be found in
[custom process](custom_process.md).

