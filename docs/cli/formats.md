# File formats

## General

By default TIPS tries to figure out the file format from its extension, with the
following rules:

- files ending with `.yml` are assumed to the PiNN tfrecord format (`-f pinn`).
- files ending with `.dump` are assumed to LAMMPS dump file (`-f dump`), they
  are handled by ASE but with some extra options (see below)
- files ending with `.trr` are assumed to be gromacs trjectoriess (`-f trr`),
  the are handelled by MDAnalysis and extra options (see below) are appectped.
- All other files are parsed to the ASE
  [`iread`](https://wiki.fysik.dtu.dk/ase/ase/io/io.html#ase.io.iread) function,
  see the ASE documentation for a list of available
  [formats](https://wiki.fysik.dtu.dk/ase/ase/io/io.html#ase.io.write). You may
  enforce the format.

Like the input function, TIPS uses ASE for file output, `-o output -of ext` will
write out try to write the dataset to `output.ext` with `ase.io.write`, with the
exeptions of `-of pinn` and `-of dump` .

## LAMMPS 

In most file formats the atoms are identified by their elements or atomic
numbers. However, the atom types are usually specified in LAMMPS as consecutive
numbers as `1, 2, 3, ...` while their elements might only be inferred from their
mass. In addition, LAMMPS has a flexible
[units](https://lammps.sandia.gov/doc/units.html) system Both make it hard to
convert a LAMMPS dump to other formats. Therefore tips provides two extra
options `units` and `emap` when reading/writing LAMMPS files.

The `emap` options specifies the mapping of elements with a string of comma
separated `in:out` pairs. For instance, suppose you have a dump file with the
elements `1:H, 2:O, 3:Na, 4:Cl` and `real` units. The following command should
convert it to an `.xyz` file with proper elements and units.

A log file in the lammps log foramt can be attached to provide the potential
energy information (it is assumed that the thermo output has the same frequency
as the dump file).

```
tips convert output.dump --emap '1:1,2:8,3:11,4:17' --units real --log output.log
```

## Gromacs

A topology file is required to read Gromacs trajectories, in addition, an energy
log file can be supplied to provide the energy labels. It is assumed that the
energy log and the trajectory file contains matching images, energy and force
labels.

```
# Assume you run the simulation like this
gmx grompp -f input.mdp -c input.pdb -p input.top -o input.tpr
gmx mdrun -s input.tpr --deffnm output
# this converts the outputs to a output.xyz file containing forces and energies
echo 4\n0 | gmx energy -f output.edr -o energy.xvg
tips convert output.trr --top input.pdb --log energy.xvg -o output -of xyz 
```
