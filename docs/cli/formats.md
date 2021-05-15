# General Information

By default TIPS tries to figure out the file format from its extension, with the
following rules:

- files end with `.yml` are assumed to the PiNN tfrecord format (`-f pinn`).
- files end with `.dump` are assumed to LAMMPS dump file (`-f dump`), they are
  handled by ASE but with some [extra options](lammps.md)
- All other files are parsed to the ASE
  [`iread`](https://wiki.fysik.dtu.dk/ase/ase/io/io.html#ase.io.iread) function,
  see the ASE documentation for a list of available
  [formats](https://wiki.fysik.dtu.dk/ase/ase/io/io.html#ase.io.write). You may
  enforce the format.

Like the input function, TIPS uses ASE for file output, `-o output -of ext` will
write out try to write the dataset to `output.ext` with `ase.io.write`, with the
exeptions of `-of pinn` and `-of dump` .

