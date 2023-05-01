#  ASE Module

The `ase.nf` module supplies two boilerplate processes `aseSP` and `aseMD` that
be used to construct processes based on ASE. A typical usage of the module is to
supply a python script that defines the an ASE calculator object as the input,
and make the SP and MD computation with ASE. Those usage can be easily achieved
by the provided `sp` and `md` workflow. The [DFTB+ module] is an example of such
usage.

[DFTB+ module]: ../dftb

## aseSP

The aseSP process perform

### Channel specification

| Element | Type   | i/o      | Note                                         |
|---------|--------|----------|----------------------------------------------|
| `name`  | `val`  | `in[0]`  | an id to identify the process                |
| `input` | `path` | `in[1]`  | a python module with a `calc` object defined |
| `geo`   | `path` | `in[2]`  | a geometry file readable by `ase.io.read`    |
| `aux`   | `path` | `in[3]`  | auxiliary file to be supplied                |
| `name`  | `val`  | `out[0]` | same as input                                |
| `sp`    | `path` | `out[1]` | single point label in xyz format [`sp.xyz`]  |
