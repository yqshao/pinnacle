# -*- coding: utf-8 -*-

from .dataset import Dataset
from .runner import ds2runner
from .pinn import ds2pinn
from .cp2k import load_cp2k
from .deepmd import load_deepmd_raw
from .lammps import load_lammps_dump


loaders = {
    'deepmd-raw': load_deepmd_raw,
    'lammps-dump': load_lammps_dump,
    'cp2k': load_cp2k
}

convertors = {
    'runner': ds2runner,
    'pinn': ds2pinn
}

def load_ds(*args, fmt='auto', **kwargs):
    if fmt in loaders:
        return loaders[fmt](*args, **kwargs)
    else:
        raise NotImplementedError(f'unknown format "{fmt}"')
