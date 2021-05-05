# Development environment 

By default TIPS processes uses singularity images. For development purpose you
might want an environment that instead follows your current code. This is can be
easily configured with the `nextflow.config` file. This document provides some
hints.

## On HPC

This section assumes that you have some modules available on the HPC you are
using, taking the [Alvis](https://www.c3se.chalmers.se/about/Alvis/) cluster as
an example.

Suppose you are working on a workflow that depends on TIPS, but you want to
tweak the TIPS code for the outputs, or you need to edit your AML code (PiNN)
here. Building an image each time you commit a change will be unnecessary, and 
also resource consuming.

Instead, we with to use the available software optimized for the cluster as much
as possible, and also a manageable development environment. Following the admin's
[suggestion](https://www.c3se.chalmers.se/documentation/applications/python/), we 
set up python a virtual environment for development.

```shell
# load available modules on HPC to avoid extra dependency
module load GCC/8.3.0  CUDA/10.1.243  OpenMPI/3.1.4\
  TensorFlow/2.3.1-Python-3.7.4 PyYAML  matplotlib
virtualenv --system-site-packages tips_env 
# install the python packages, suppose we work on both TIPS and PiNN
git clone https://github.com/yqshao/tips.git
git clone https://github.com/yqshao/PiNN.git -b TF2
pip install -e tips/python
pip install -e PiNN
```

Then, load the modules in the `nextflow config`

```groovy
profiles {
  default {
    process {
      executor = 'slurm'
      params.lmpCmd = 'singularity exec docker://lammps/lammps lmp_serial'
      module = 'GCC/8.3.0:CUDA/10.1.243:OpenMPI/3.1.4:TensorFlow/2.3.1-Python-3.7.4:PyYAML:matplotlib:GROMACS'
      beforeScript = 'source $HOME/tips_env/bin/activate'
      clusterOptions = '--gres=gpu:T4:1'
      withLabel: tips {
        executor = 'local'
      }
      withLabel: pinn {
        time='12h'
      }
      withLabel: lammps {
        time='2h'
      }
      withLabel: gromacs {
        time='2h'
      }
    }
  }
}
```

- The moduels will be loaded and you virtual environment enabled for you processes
- Extra binaries can still be specified (e.g. with a `lmpCmd` here)
- You can fine tune the queuing option for each process

Now you will be able to `nf run your_script.nf`, and your own version of PiNN
and TIPS will be used (without images built).
