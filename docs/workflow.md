# Workflows

This document introduces workflows implemented in TIPS, long with templates for
those workflows.

## Continuous training and sampling

In a typical learning and sampling scheme, the model is continuously evaluated
by new samples, and models are trained until the model is considered as
converged for the newly sampled data.

!!! note 
    Here and in most cases, a dataset can be replaced with a sampler and a sampler
    with a `Channel` of datasets (e.g. by splitting a dataset). TIPS will try to 
    generate the `initiDs` with the sampler if one is not given.

``` mermaid
graph LR
  A[initDs] --> B([trainer]);
  B --> C[model]:::inter;
  D([sampler]) --> E{{tol?}};
  C --> E;
  E -- yes --> F[finalModel];
  E -- no --> G[augDs];
  G --> B;
```

### Example

### Options 


## Active Learning

In an active learning workflow, the model is used as a input for the sampler.

``` mermaid
graph LR
  A[initDs] --> B([trainer]);
  B --> C[model];
  C --> D([sampler]);
  D --> E{{tol?}} 
  E -- yes --> F[finalModel];
  E -- no --> G[augDs];
  G --> B;
```

### Example

### Options 
