# -*- coding: utf-8 -*-

from .dataset import Dataset
from .runner import ds2runner
from .pinn import ds2pinn
from .deepmd import load_deepmd_raw


loaders = {
    'deepmd-raw': load_deepmd_raw
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
