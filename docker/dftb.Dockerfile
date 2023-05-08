FROM mambaorg/micromamba:1.4.2-bullseye
RUN micromamba install --yes --name base --channel conda-forge \
    numpy==1.22.1 ase==3.22.1 dftbplus=22.2=mpi_openmpi_* \
    && micromamba clean --yes --all
ENV PATH "$MAMBA_ROOT_PREFIX/bin:$PATH"
CMD python
