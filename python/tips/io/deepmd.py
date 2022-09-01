# -*- coding: utf-8 -*-

"""The DeepMD dataet loader

Originally written by: Alexandros Zantis
Adapted by: Yunqi Shao [yunqi.shao@kemi.uu.se]
"""
import numpy as np
from tips.io.utils import list_loader


@list_loader
def load_deepmd_raw(directory, eunit=1.):
    """The raw format specification:
    https://docs.deepmodeling.com/projects/deepmd/en/latest/data/data-conv.html#raw-format-and-data-conversion
    Note the at DeePME-kit does not fix a standard unit for energy, and
    the energy will be read as is, if not specified.

    Args:
        directory (str or list of str): folder for the raw files
        eunit: energy unit

    Returns:
        Dataset: a TIPS dataset
    """
    from ase.data import atomic_numbers
    from tips.io.dataset import Dataset

    # the naming scheme is the default one
    energy = np.loadtxt(f"{directory}/energy.raw")*eunit
    nlabel = energy.shape[0]
    force = np.loadtxt(f"{directory}/force.raw").reshape([nlabel, -1, 3])*eunit
    coord = np.loadtxt(f"{directory}/coord.raw").reshape([nlabel, -1, 3])
    cell = np.loadtxt(f"{directory}/box.raw").reshape([nlabel, 3, 3])
    type = np.loadtxt(f"{directory}/type.raw", dtype=np.int)
    type_map = np.genfromtxt(f"{directory}/type_map.raw", dtype="str")

    # cast the dataset to a universal form
    type_map = np.array([atomic_numbers[type] for type in type_map], dtype=np.int)
    elem = np.tile(type_map[type], [nlabel, 1])

    all_data = {
        "cell": cell,
        "elem": elem,
        "coord": coord,
        "force": force,
        "energy": energy,
    }

    # generating meta data
    meta = {
        "fmt": "DeepMD raw",
        "size": coord.shape[0],
        "elem": set(type_map),
        "spec": {
            "cell": {"shape": [3, 3], "dtype": "float"},
            "elem": {"shape": [None], "dtype": "int"},
            "coord": {"shape": [None, 3], "dtype": "float"},
            "force": {"shape": [None, 3], "dtype": "float"},
            "energy": {"shape": [], "dtype": "float"},
        },
    }

    def indexer(i):
        datum = {k: v[i] for k, v in all_data.items()}
        return datum

    return Dataset(meta=meta, indexer=indexer)
