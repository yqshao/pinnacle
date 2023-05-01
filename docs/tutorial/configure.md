# Profiles and Configuration

## The nextflow.config file

In your workflow folder, you will find a file named `nextflow.config` that
looks like this:

```groovy
--8<-- "nextflow.config::3"
```

The `nextflow.config` file contains detailed specifications about the resources
to be used by each process, it can also change default parameters in the
workflow. The configuration file separates the definition from that of the
workflow itself, which means the same workflow will be runnable on any resources
supported by nextflow (from local computer, to HPC computer, to cloud-based
platforms). Available options can be found in the nextflow [documentation].

[documentation]: https://www.nextflow.io/docs/latest/config.html

In PiNNAcLe, those configurations are always defined as "profiles", so that they
can be switched with the `-profile name` option for `nextflow run`. PiNNAcLe
curates a list of [profiles](../../profiles/overview) for some of the
computational resources we have access to. Those can be good starting points for
you to adapt for your own need. You are also welcome to contribute your config
if think it will be useful for others.
