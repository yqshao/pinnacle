# Reading LAMMPS files

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

```
tips convert input.dump --emap '1:1,2:8,3:11,4:17' --units real
```

