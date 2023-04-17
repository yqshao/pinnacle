# Alvis

> The [Alvis] cluster is a national NAISS resource dedicated for Artificial
> Intelligence and Machine Learning research.

[Alvis]: https://www.c3se.chalmers.se/about/Alvis/


## TensorFlow and PiNN

Since TensorFlow is already installed on Alvis, it's recommended to run PiNN
with that. To do so, make a python environment with the supplied TF

```bash
ml TensorFlow/2.5.0-fosscuda-2020b
python -m venv $HOME/pinn-tf25
source $HOME/pinn-tf25/bin/activate
pip install git+https://github.com/teoroo-cmc/pinn.git 
```

The above creates a python virtual environment based on the system TF module.
When you need to run PiNN manually in a new bash session, you need to load the
module and activate the environment:

```
ml TensorFlow/2.5.0-fosscuda-2020b
source $HOME/pinn-tf25/bin/activate
```

## CP2K

The container image in [NGC] for CP2K supports acceleration through CUDA. You
will need to build the singularity file following the NGC instructions, and
point the profile to your image. The accelerators should be picked up
automatically, for which you can verify by looking for the `ACC:` tags in the
CP2K log file.

[NGC]: https://catalog.ngc.nvidia.com/orgs/hpc/containers/cp2k

## Profile

```groovy
  alvis {
    executor.name = 'slurm'
    executor.queueSize = 100
    executor.submitRateLimit = '120 min'
    params.cp2k_cmd = 'OMP_NUM_THREADS=2 mpirun -n 4 cp2k.psmp'
    process {
      time = '3d'
      executor = 'slurm'
      errorStrategy='ignore'
      beforeScript = 'source $HOME/pinn-tf25/bin/activate'
      module = 'TensorFlow/2.5.0-fosscuda-2020b'
      withLabel: 'tips|moltemplate' {executor='local'}
      withLabel: latent {
        time = '1h'
        clusterOptions = '--gres=gpu:T4:1 -A YOUR_PROJ_ID'
      }
      withLabel: 'pinn|ase' {
        scratch=true
        clusterOptions = '--gres=gpu:T4:1 -A YOUR_PROJ_ID'
      }
      withLabel: 'cp2k' {
        scratch=true
        clusterOptions = '--gres=gpu:T4:2 -A YOUR_PROJ_ID'
        container='$HOME/cp2k_v9.1.0.sif'
      }
    }
  }

```
