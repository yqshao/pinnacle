# The TIPS IO module

The IO module in TIPS allows for the translation between different atomistic
data formats, with a special focus for AML. The core class is `tips.io.Dataset`,
which holds the dataset and its metadata, and allows for manipulation and
convertion of the data. To create datasets, the easiest way is to use the
universal data loader `load_ds`:

=== "Usage"

    ```Python
    from tips.io import load_ds
    
    ds = load_ds('cp2k-pos-1.xyz',
                 frc='cp2k-frc-1.xyz',
                 cell='cp2k-cell.dat',
                 fmt='cp2k')

    # data can be splitted, converted
    ds.split({'train':8, 'test':2})
    ds.convert('export.xyz', fmt='ext-xyz')
    print(ds)
    ```


=== "Output"

    ```Python
    # printing the  dataset shows basic information about the dataset
    <tips.io.Dataset:
     size: 50
     format: 'CP2K outputs'
     indexable: True
     structure:
       - elems: ?, int32
       - coord: ?x3, float32
       - cell: 3x3, float32>
    ```

Detailed descriptions about the `Dataset` object can be found in its [API documentation]().


## Available formats

!!! warning "Check marks are planned implementations for now"

| Format  | Read               | Convert            | Note                           |
|---------|--------------------|--------------------|--------------------------------|
| ase     | :white_check_mark: | :white_check_mark: | ASE Atoms objects              |
| cp2k    | :white_check_mark: |                    | CP2K data (pos, frc, and cell) |
| deepmd  | :white_check_mark: | :white_check_mark: | DeePMD format                  |
| ext-xyz | :white_check_mark: | :white_check_mark: | Extended XYZ format            |
| lammps  | :white_check_mark: | :white_check_mark: | LAMMPS dump format             |
| runner  | :white_check_mark: | :white_check_mark: | RuNNer format                  |

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

## Registered datasets

TIPS curates a small list of datasets that can be directly accessed via the
`load_ds` function. For now, the list exist mainly for test and demonstrative
purpose.

!!! warning "Not Implemented yet!"

## API documentation

!!! warning "Not Implemented yet!"
