# Understanding the workflow

## Parameters

So far in the tutorial, we have executed the workflow as-is. In nextflow, there
is a built-in way of defining variables that can be tweaked at runtime, called
"parameters". See the `nextflow/we_demo.nf` as an example:

```groovy
--8<-- "nextflow/we_demo.nf:params"
```

The `params.proj` parameter determines where the output files will be
"published" into. To change while running the same workflow, we simply have to
add the option to `nextflow run`, i.e. the following command will run the same
workflow, but use 2 initial seeds (to generate the initial dataset with).

```bash
nextflow run main.nf --proj=trail2 --init_seeds=2
```

## Modules and inclusion

Thanks to the DSL2 syntax of nextflow, the workflows in PiNNAcLe are reusable
easily. Below, on see how the workflow parameters can be injected into included
processes/workflows; the publishing directories of subworkflows is set, and a
few parameters required by acle are supplied:

```groovy
--8<-- "nextflow/we_demo.nf:include"
```


## Processes 

The processes, such as the ones included above are the building block of
nextflow scripts. As shown below, they are simply scripts with a few
input/output "channels". When the workflow runs, the scripts will be executed in
the work directory, where the `$variables` will be substituted by whatever
supplied by the input channels; and the files produced by the output patterns
will.

```groovy
--8<-- "nextflow/module/ase.nf:sp"
```

It should be easy to tweak the workflow behavior via the process scripts. You
might also want to check the nextflow [documentation] to see what is possible.
In addition, a few implementation patterns are followed:

[documentation]: https://www.nextflow.io/docs/latest/process.html

- processes are labelled by the dependent software and arranged as "modules";
- interchangeable processes in different modules share the same name (e.g. `sp`,
  `md`);[^1]
- processes input is always one tuple, where the first element is indented as an
  identifier;[^2]
  
## Workflow and recipes

IN DSL2, workflows are composed of processes, and can share the same
[input/output][wf_io] structure to behave like a process. However, the "main"
workflows to be executed would have to be self-contained, meaning it will not
receive inputs and will define the "channels" within itself and feed them to
the sub-workflows, as shown in the following snippet.

[wf_io]: https://www.nextflow.io/docs/latest/dsl2.html#workflow-input

```groovy
--8<-- "nextflow/we_demo.nf:wf"
```

Those modules that provide self-contained workflow are called "recipes" in
PiNNAcLe. The workflow will always be named as `entry` and exposed as the module
name in the main workflow, i.e., the following three commands does the same
thing:

```bash
nextflow run nextflow/acle.nf
nextflow run nextflow/acle.nf -entry entry
nextflow run main.nf -entry acle
```

I.e., the `-entry` option chooses which workflow to run within a workflow file,
and the `main.nf` serves as a shortcut to call those recipes. Other than that,
there are no qualitative difference between recipes and normal modules, and both
can expose reusable sub-workflows, such as the `acle` workflow in the demo.


[^1]: This also applies to workflows that shares a common input/output pattern.
    This allows complex workflows to be composed in a implementation-agnostic
    fashion, for instance, the reference (`params.ref`) and machine learning
    potential (`params.mlp`) can be swapped to any compatible module.

[^2]: This ensures all outputs are trackable by their identifier in a complex
    workflow. For instance, in the acle workflow, the labelled points are
    identified as `gen/geo/idx`, the information is used to map them back to the
    sampled trajectory when evaluating the performance of the model.
