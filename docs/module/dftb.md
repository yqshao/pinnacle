#  DFTB Module

The `dftb.nf` module is an alias to the ASE module that provides the `sp` and
`md` process. An example input file for both processes (`input/dftb/xtb.py`) is
shown below:

```python
--8<-- "input/dftb/xtb.py"
```


The label and aux files are redirected to `'dftb'` and `params.dftb_aux`.
See the [ASE module] documentation for more details. 

[ASE module]: ../ase
