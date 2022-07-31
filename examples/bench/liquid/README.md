# Example benchmark workflow for aqueous systems

This is an example workflow with TIPS for aqeuos systems, this workflow trains a
set of MLP models and test them with NPT simulations. The workflow does not
contain any dataset, to run it, download one of the datasets and (and convert it
to a recognizable format if necessary). To run this example, consider copying
the required files to your working directory

## Prepare your environment

This workflow is tested on the Alvis cluster and the default config is composed
for Alvis, to set up a development environment:

```bash
ml Nextflow/20.10.0 TensorFlow/2.5.0-fosscuda-2020b
virtualenv --system-site-packages ~/tips-env
pip install git+https://github.com/Teoroo-CMC/PiNN.git
git clone git@github.com:yqshao/tips.git
pip install -e tips
```

Run the following before running the scripts:

```bash
ml Nextflow/20.10.0 TensorFlow/2.5.0-fosscuda-2020b
source ~/tips-env/bin/activate
```

To run this example, consider copying the files to a new folder

```bash
cp -rL $TIPS_INSTALL_DIR/examples/bench/liquid my_liquid_bench
cd my_liquid_bench
```

Optionally, you can copy one of our HPC profiles as a template for your local
HPC systems, e.g.:

```bash
cp -rL $TIPS_INSTALL_DIR/profiles/alvis.config ./nextflow.config
```

## Getting some data

To get the [HDNNP](https://zenodo.org/record/2634098) training data for H2O:

```bash
wget -qO- https://zenodo.org/record/2634098/files/training-data_H2O.tar.gz?download=1 | tar xvz
```

For a larger
[dataset](https://figshare.com/articles/dataset/Data_from_Dissolving_salt_is_not_equivalent_to_applying_a_pressure_on_water_/17193023/1)
from [Zhang et al., 2021](https://www.nature.com/articles/s41467-022-28538-8):

```bash
wget -q https://figshare.com/ndownloader/files/31775396 && unzip 31775396 && rm 31775396
tips convert data/taining_data/*/* -f deepmd-raw -o salt.data -of runner
```

## Run the benchmark

To run the benchmark:

```bash
nextflow run main.nf # optionally: -profile alvis
```

## What's next

Now that you have a workflow running, there are several things you can do to
extend this workflow, you might want to check how to:

- [including custom analysis]()
- [including new MLP models]()
- [tweaking the workflow]()
