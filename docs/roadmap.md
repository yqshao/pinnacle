# Roadmap

TIPS is developed as a open-source project and is planned to be released under
the BSD 3-clause license. However, the code is at an early stage of development,
**please use at your own risk**. 

Below are some known limitations of TIPS in case you are going to proceed. There
will be an official release after we solve those issues (or decide to live with
them).
    
## Known limitations

- TIPS relay on a hackish
   [implementation](https://github.com/yqshao/tips/blob/c9f9f6b83ce19342074089aff3fd9e0f876b5cb8/nextflow/utils.nf#L18-L34)
   of feedback loops (see also [this
   issue](https://github.com/nextflow-io/nextflow/issues/1766)).
- The way TIPS handles inputs partially breaks the caching mechanism of Nextflow.
- Too little examples (more examples on reusing the workflows are planned).
- File formats, currently all datasets are converted to `.xyz` files between processes, it works Ok for potentials, but we might want to think about if we want to use other (more generic or flexible) formats.

Some of those issues might be solved by extending the Nextflow Channel and
process classes.
