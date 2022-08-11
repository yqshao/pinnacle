# Utility processes

## convert

The convert process converts one dataset to another, where the `dataset` input
can be a list of files (in such case), the dataset will be joined together.

| Channel      | Type | Note                          |
| ------------ | ---- | ----------------------------- |
| (in) name    | val  | an id to identify the process |
| (in) dataset | file | input dataset(s)              |
| (in) flags   | val  | `tips convert` options        |
| (out) output | file | the converted dataset         |
