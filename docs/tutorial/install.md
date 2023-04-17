# Get started

## Install TIPS via pip

TIPS is a chain of tools towards the construction of machine-learnt interatomic
potentials. At its core, TIPS consists of a hierachy of [nextflow]() workflows
for a variety of atomistic machine learning tasks. If you are familiar with
nextflow, you can run the [recipes]() without even installing anything. For
newcomers, it's recommanded to install the TIPS CLI, via the `pip` command:

=== "Command"

    ```bash
    pip install 'git+https://github.com/yqshao/tips.git#egg=tips&subdirectory=python'
    ```

## First step

To get started, the easiest way is to create a new project with the
interactive wizard command:

=== "Command"

    ```bash
    tips wizard
    ```

The wizard will walk you through the configuration of your workflow and create a
folder with the necessary files. To run the worflow, you'll also need to install
Nextflow (the wizard will direct you to the nextflow installation guide if it is
not found). The next sections will detail the contents in the folder and the
execution of the workflow.
