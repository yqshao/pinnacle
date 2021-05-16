# Command Line Interface (CLI)

TIPS provides a command interface as gluing utilities between different trainer
and sampler modules.

## Dependencies

- click >= 7.0
- numpy >= 1.3.0
- ase >= 3.19.0
- pyyaml >= 3.01

In addition, tensorflow >= 2.0 is required to read and write tfrecord files used
by PiNN and MDAnalysis required for reading Gromacs trr files.

## Installation

The CLI of TIPS is implemented as a python package, to install it:

```
pip install git+https://github.com/yqshao/tips#subdirectory=python
```

See the documentation or run `tips -h` to see a list of availabel commands.
