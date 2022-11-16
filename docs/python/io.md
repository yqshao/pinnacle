# The TIPS IO module

The IO module in TIPS allows for the translation between different atomistic
data formats, with a special focus for AML. 

## Available formats

| Format      | Read               | Convert            | Note                           |
|-------------|--------------------|--------------------|--------------------------------|
| ase/asetraj | :white_check_mark: | :white_check_mark: | ASE Trajectory obj or files    |
| cp2k        | :white_check_mark: | :no_entry_sign:    | CP2K data (pos, frc, and cell) |
| deepmd      | :white_check_mark: | :no_entry_sign:    | DeePMD format                  |
| extxyz      | :no_entry_sign:    | :white_check_mark: | Extended XYZ format            |
| lammps      | :white_check_mark: | :no_entry_sign:    | LAMMPS dump format             |
| pinn        | :white_check_mark: | :white_check_mark: | PiNN-style TFRecord format     |
| runner      | :no_entry_sign:    | :white_check_mark: | RuNNer format                  |


## General specification

### units and formats

TIPS uses a unit system compatible to ASE internally, that is:

- energy in eV
- length in Å

Some formats does not have a fixed unit system, or a different unit standard,
those are documented in the format-specific informaiton section

### `load_ds` funciton

The `load_ds` function is a universal entry point for dataset loaders in TIPS.

=== "Usage"

    ```Python
    from tips.io import load_ds

    ds = load_ds('path/to/dataset', fmt=deepmd-raw')
    ds = ds.join(ds) # datasets can be joined together
    ds.convert('dataset.yml', fmt='pinn') # the and converted to different formats
    print(ds)
    ```

=== "Output"

    ```Python
    # printing the  dataset shows basic information about the dataset
    <tips.io.Dataset:
     fmt: DeePMD raw
     size: 100
     elem: 8, 1
     spec:
         cell: 3x3, float
         elem: [None], int
         force: [None, 3], float
         coord: [None, 3], float
         energy: [], float>
    ```

The function returns a `Dataset`-class object, its usage is detailed below.

### `Dataset` class
::: tips.io.dataset.Dataset


## Format specific info.

### tips.io.ase

The `tips.io.ase` module allows the conversion of ase trajectories into the TIPS
Dataset. Writer for `asetraj` and `extxyz` extends the original ASE file
writers, and adds additional columns such as `force_std` which one might obtain
in an ensemble-based MD simulation.

??? "Source code"

    ```python
    --8<-- "python/tips/io/ase.py"
    ```

### tips.io.cp2k

The `cp2k` module reads CP2K ouputs in written as specified in the
[`%MOTION%PRINT%`](https://manual.cp2k.org/trunk/CP2K_INPUT/MOTION/PRINT.html)
section. Those files are typically named as `path/proj-pos-1.xyz`,
`path/proj-frc-1.xyz`, etc, where the project name are specified in
[`%GLOBAL%PROJECT_NAME%`](https://manual.cp2k.org/trunk/CP2K_INPUT/GLOBAL.html#list_PROJECT_NAME).
**Note** that the CP2K format assumes atomic units, and loader uses CODATA
version 2006, as adapted by CP2K instead of the 2014 version used in ASE by
default.

??? "Source code"

    ```python
    --8<-- "python/tips/io/cp2k.py"
    ```

### tips.io.cp2klog

The `cp2klog` module reads in information as specified in
[`%FORCE_EVAL%PRINT%`](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/PRINT.html)
section. Those outputs will by wriiten by CP2K to stdout. Like the `cp2k`
module, `cp2klog` assumes the atomic units and uses CODATA 2006.

??? "Source code"

    ```python
    --8<-- "python/tips/io/cp2klog.py"
    ```

### tips.io.lammps

The `lammps` reads the lammps formatted `.dump` files, note that this
implementation only supports the limited format with the atom format: `ITEM:
ATOMS id type x y z`, any other format should fail with an error. For lammps
files it's common that the "real" elements information is stored in a separate
`.data` file. The element can be converted easily with the `.map_elems()` method
of the `Dataset` class.

??? "Source code"

    ```python
    --8<-- "python/tips/io/lammps.py"
    ```

### tips.io.runner

The runner file format uses atomic units.

??? "Source code"

    ```python
    --8<-- "python/tips/io/runner.py"
    ```

## Custom reader/writer

It is possible to extend TIPS by registering extra reader/writers, an example
for custom reader/writer can be found below:

??? example "Example implementation of custom data reader and converters"

    ```Python
    from tips.io.utils import tips_reader, tips_convert

    @tips_reader('my-ase')
    def load_ase(traj):
        """ An example reader for ASE Atoms

        The function should return a tips Dataset, by specifying at least an generator which
        yields elements in the dataset one by one, and the metadata specifying the data structure.
        (the generator is redundent ni the below case because an ASE trajectory is indexable
        and has a defined size, such a generator will be defined automatically by tips)

        Args:
            traj: list of atoms

        Returns:
            tips.io.Dataset
        """
        from tips.io import Dataset
        meta = {
          'spec': {
            'elems': {'shape': [None], 'dtype': 'int32'},
            'coord': {'shape': [None, 3], 'dtype': 'float32'},
            'cell': {'shape': [3, 3], 'dtype': 'float32'}
          },
          'size': len(traj),
          'fmt': 'Custom ASE Format'
        }

        def indexer(i):
            atoms = traj[i]
            data = {
               'elems': atoms.numbers
               'coord': atoms.positions,
               'cell': atoms.cell,
            }
            return data

        def generator():
            for i in range(meta['size']):
                yield indexer(i)

        return Dataset(generator=generator, meta=meta, indexer=indexer)


    @tips_convert('my-ase')
    def ds_to_ase(dataset):
        """ An example data converter to ASE trajectory

        The function must takes on dataset and optionally extra keyword arguments as inputs.
        There is no limitaton on the return values.

        Args:
            dataset (tips.io.Dataset): a dataset object
        """
        from ase import Atoms
        traj = [
           Atoms(data['elems'],
                 positions=data['coord'],
                 cell=data['cell'])
           for data in dataset
        ]
        return  traj
    ```

    The additonal format will be available for data loading and conversion:
    ```Python
    ds = load_ds([Atoms['H'], Atoms['Cu']], fmt='my-ase')
    traj = ds.convert(fmt='my-ase')
    ```

