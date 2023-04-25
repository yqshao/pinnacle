# The Alvis Cluster

The [Alvis] cluster is a national NAISS resource dedicated for Artificial
Intelligence and Machine Learning research. The Alvis profile is configured to
use the GPU resources there.

[Alvis]: https://www.c3se.chalmers.se/about/Alvis/

## TensorFlow and PiNN

Since TensorFlow is already installed on Alvis, it's recommended to run PiNN
with that. To do so, make a python environment with the supplied TF:

```bash
ml TensorFlow/2.6.0-foss-2021a-CUDA-11.3.1
python -m venv $HOME/pinn-tf26
source $HOME/pinn-tf26/bin/activate
pip install git+https://github.com/teoroo-cmc/pinn.git 
```

The above creates a python virtual environment based on the system TF module.
When you need to run PiNN manually in a new bash session, you need to load the
module and activate the environment:

```bash
ml TensorFlow/2.6.0-foss-2021a-CUDA-11.3.1
source $HOME/pinn-tf26/bin/activate
```

You might also want to make this enivronment avaialble to the [Alvis
OnDemand][ondemand] portal, following the instruction (after activating your
environment):

```bash
pip install ipykernel
python -m ipykernel install --user --name pinn-tf26 --display-name "pinn-tf26"
```

[ondemand]: https://portal.c3se.chalmers.se/pun/sys/dashboard/

## CP2K

The container image in [NGC] for CP2K supports acceleration through CUDA. You
will need to build the singularity file following the NGC instructions, and
point the profile to your image. The accelerators should be picked up
automatically, for which you can verify by looking for the `ACC:` tags in the
CP2K log file.

[NGC]: https://catalog.ngc.nvidia.com/orgs/hpc/containers/cp2k

## Profile

```groovy
--8<-- "profiles/alvis.config"
```
