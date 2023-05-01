# Profiles for PiNNAcLe Workflows

This part of the documentation hosts the setups to run PiNNAcLe workflows on
different computational resources. The instructions also include the
installation of libraries to run optimally on specific resource.

!!! warning ""

    Please note that the softwares installed on supercomputer are subject to changes, 
    and the developers are not tracking all of them. You are nevertheless welcome to 
    report deprecated instructions or update them.

## Standard profile

```groovy
--8<-- "profiles/standard.config"
```

The standard profiles uses container images for its processes, most of them are
compiled hosted on Dockerhub and automatically built (with a few exceptions).
Below are a list of available softwares and links to the docker image, where the
build scripts are also kept.

| Software           | Version         | Dockerhub image                                      |
|--------------------|-----------------|------------------------------------------------------|
| packmol            | 20.14.2         | [teoroo/pinnacle:molutils]                           |
| molptemlate        | 2.20.19         | [teoroo/pinnacle:molutils]                           |
| awk, bc, openbabel | latest (Debian) | [teoroo/pinnacle:molutils]                           |
| dftbplus           | 22.2            | [teoroo/pinnacle:dftbplus]                           |
| tips               | latest          | [teoroo/pinnacle:pinn]                               |
| pinn               | latest          | [teoroo/pinnacle:pinn]                               |
| cp2k               | 9.1.0 (NGC)     | [nvcr.io/hpc/cp2k:v9.1.0]                            |
| lammps             | 29Sep2021       | [lammps/lammps:stable_29Sep2021_centos7_openmpi_py3] |

[teoroo/pinnacle:molutils]: https://hub.docker.com/r/teoroo/pinnacle/molutils
[teoroo/pinnacle:dftbplus]: https://hub.docker.com/r/teoroo/pinnacle/dftbplus
[teoroo/pinnacle:tips]: https://hub.docker.com/r/teoroo/pinnacle/tips
[teoroo/pinnacle:pinn]: https://hub.docker.com/r/teoroo/pinnacle/pinn
[nvcr.io/hpc/cp2k:v9.1.0]: https://catalog.ngc.nvidia.com/orgs/hpc/containers/cp2k/tags
[lammps/lammps:stable_29Sep2021_centos7_openmpi_py3]: https://hub.docker.com/r/lammps/lammps/stable_29Sep2021_centos7_openmpi_py3
