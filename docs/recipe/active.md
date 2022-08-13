# Active workflows

TIPS implements several strategies to run active learning workflows, as detailed
below.

## Query by Committee

The query by committee (QbC) takes one dataset and adaptively samples a small
fraction of data from it. Note that unlike the other AL workflows, this scheme
does not actively generate new dataset.

=== "Flowchart"

    ```mermaid
    graph LR
    filter([Filter]) --> ds[Dataset]
    ds ------ |End or Next Iter.| qbc
    ref[Reference] --> filter
    ref -----> qbc([QbC])
    inp[Input] ---> train([Train])
    ds --> train
    seeds[Seeds] ---> train
    subgraph QbC Iteration
      train --> model[Model]
      model --> qbc
    end
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

## Ensemble NN

In an ensemble NN workflow, the trajectory is propagated with an ensemble of NN
models, while a subset of the trajectory is labelled according to a given
uncertainty tolerance.

=== "Flowchart"

    ```mermaid
    graph LR
    ds[Dataset] --- |End or Next iter.| filter([Filter])
    ds & inp[Input] & seeds[Seeds] --> train([Train])
    subgraph ENN Iteration
      train --> model[Model]
      model --> emd([Ensemble MD])
      emd ---> traj[Trajecotry] & Uncertainty --> filter
    end
    ```

## Density-based Clustering

The density-based clustering (DBC) scheme samples the configuration space by
actively biasing the dynamics according to the data distribution in latent
space.

## Other links

- [Sample model input files]()
- [Adding a custom analysis]()
