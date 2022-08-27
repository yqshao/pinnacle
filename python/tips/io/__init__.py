# -*- coding: utf-8 -*-

from .dataset import Dataset
from .ase import load_ase, load_asetraj, ds2ase, ds2asetraj, ds2extxyz
from .cp2k import load_cp2k
from .cp2klog import load_cp2klog
from .deepmd import load_deepmd_raw
from .lammps import load_lammps_dump
from .pinn import ds2pinn
from .runner import ds2runner

loaders = {
    'ase': load_ase,
    'asetraj': load_asetraj,
    'cp2k': load_cp2k,
    'cp2klog': load_cp2klog,
    'deepmd-raw': load_deepmd_raw,
    'lammps-dump': load_lammps_dump,
}

convertors = {
    'ase': ds2ase,
    'asetraj': ds2asetraj,
    'extxyz': ds2extxyz,
    'pinn': ds2pinn,
    'runner': ds2runner,
}

def load_ds(*args, fmt='auto', **kwargs):
    if fmt in loaders:
        return loaders[fmt](*args, **kwargs)
    else:
        raise NotImplementedError(f'unknown format "{fmt}"')
