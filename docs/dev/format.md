# Adding a dataset format

To add a file format, you need to extend the `DatasetReader` and `DatasetWriter`
class in the `tips.io` modules. The classes define the availability of labels
and how they can be read/written.

## The DataReader class

A DataReader should at least implement a `__init__` method that defines data
formats given a input `dataset`, optionally, if the format can have variable
units, a unit argument can be supplied:

```python
class MyDs(DatasetReader):
    def __init__(self, dataset, unit=None):
        # comment attributes in `DataReader`s
        self.indexable = True
        self.size = None
        self.ds_spec = {}
        # reserved for `.iter()` or `.index()` methods
        self._dsfile = dataset
        self._index = None
```

The `indexable` and `size` attribute determines how the dataset might be read,
when `indexable` is `False` and size is not available, the dataset will only be
read iteratively, this makes tasks like splitting the dataset more expansive,
since the dataset must be enumerated once before the split can be determined.

Whenever possible, it is advisable to implement the `.index()` function if a
specific data point can be retrieved without loading the entire dataset. For
plain text formats, this is possible by a fast scan through the file for certain
pattern and subsequently rewind the text file with the `file.seek()` function,
an example is the `runner.RunnerLoader`.

## The DataWriter class

A DataWriter is also initialized with the `dataset` specifying the output
location and the optional `unit` arugment. The `.write()` method should write
one sample to the dataset, and the `.finalize()` method to finalize the writing.
When implementing the `write()` method, it is again advisable to avoid caching
the entire dataset in the memory.

## Reserved labels

Some generally used label names are hard-coded, which have fixed dimension
definitions.

| Name     | Shape         | Type    | Description                    |
| -------- | ------------- | ------- | ------------------------------ |
| `coord`  | `[natoms, 3]` | `float` | cartesian coordinates of atoms |
| `force`  | `[natoms, 3]` | `float` | atomic forces                  |
| `charge` | `[natoms, 3]` | `float` | atomic chages                  |
| `elem`   | `[natoms]`    | `float` | atomic numbers                 |
| `cell `  | `[3, 3]`      | `float` | cell vectors                   |
| `pbc `   | `[3]`         | `bool`  | periodic boundary condition    |
| `energy` | `[]`          | `float` | total energy                   |

When those data exist in the dataset, they should be named consistently.

## Checklist

Below is a checklist for adding a new file format to the TIPS repo:

- [ ] Implement the format in `python/io/your_format.py`;
- [ ] Document the file format in `docs/python/io.md`;
- [ ] Adding a unit test case in `python/test/io.py`;
- [ ] Submit your pull request!

To make a smooth PR, also kindly check the developer setup and contributing
guideline.
