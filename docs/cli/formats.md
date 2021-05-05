# General Information

## Inputs

By default TIPS tries to figure out the file format from its extension, with the
following rules

- files end with `.yml` are assumed to PiNN formats.
- files end with `.dump` are assumed to LAMMPS dump file.
- All other files are handled by the ASE [`iread`](https://wiki.fysik.dtu.dk/ase/ase/io/io.html#ase.io.iread) function

## Outputs

At this point TIPS supports three types of output formats:

- PiNN tfrecrds (`-of pinn`)
- LAMMPS dump files (`-of lammps`)
- extended-xyz files (`-of xyz`)

