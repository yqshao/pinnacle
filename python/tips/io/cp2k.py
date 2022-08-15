# -*- coding: utf-8 -*-

"""The CP2K data loader"""
import numpy as np
import logging
from tips.io.utils import list_loader

logger = logging.getLogger('tips')

def _index_xyz(fname):
    """Indexes a list of cp2k files, return the location to each frame"""
    import mmap, re

    f = open(fname, "r")  # use the first line as the frame identifier
    first_line = f.readline()
    f.seek(0)
    regex = str.encode("(^|\n)" + first_line[:-1] + "(\r\n|\n)")
    m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
    locs = [match.span()[-1] for match in re.finditer(regex, m)]
    indexes = list(zip([fname] * len(locs), locs))
    f.close()
    return indexes


def _index_cell(fname):
    """Reads a cp2k cell file, returns cell vectors for all frames"""
    cells = np.loadtxt(fname, usecols=[2,3,4,5,6,7,8,9,10]).reshape([-1,3,3])
    return cells


def _index_ener(fname):
    """Indexes a cp2k energy file, returns the total energy"""
    energies = np.loadtxt(fname, usecols=4)
    return energies


def _load_pos(index):
    from ase.data import atomic_numbers

    fname, loc = index
    elem = []
    coord = []
    f = open(fname, "r")
    f.seek(loc)
    f.readline()
    while True:
        line = f.readline().split()
        if len(line) <= 1:
            break
        elem.append(atomic_numbers[line[0]])
        coord.append(line[1:4])
    f.close()
    return {"elem": np.array(elem, np.int), "coord": np.array(coord, np.float)}


def _load_frc(index):
    from ase.data import atomic_numbers

    fname, loc = index
    force = []
    f = open(fname, "r")
    f.seek(loc)
    f.readline()
    while True:
        line = f.readline().split()
        if len(line) <= 1:
            break
        force.append(line[1:4])
    f.close()
    return {"force": np.array(force, np.float)}


def _load_ener(energy):
    return {"energy": energy}


def _load_cell(cell):
    return {"cell": cell}


_frc_spec = {"force": {"dtype": "float", "shape": [None, 3]}}
_cell_spec = {"cell": {"dtype": "float", "shape": [3, 3]}}
_ener_spec = {"energy": {"dtype": "float", "shape": []}}
_pos_spec = {
    "elems": {"dtype": "int", "shape": [None]},
    "coord": {"dtype": "float", "shape": [None, 3]},
}


@list_loader
def load_cp2k(project,
              cp2k_pos="pos-1.xyz", cp2k_cell="cell-1.cell",
              cp2k_frc="frc-1.xyz", cp2k_ener="1.ener"):
    """Loads cp2k-formatted outputs as datasets

    By default, the following supported output files are scanned and matching
    files will be loaded. The loader assumes that all files contains the same
    number of matching frames.

    - ener: f"{project}-{cp2k_ener}"
    - cell: f"{project}-{cp2k_cell}"
    - pos: f"{project}-{cp2k_pos}"
    - frc: f"{project}-{cp2k_frc}"

    TODO (extra keywords) to be scanned:
    - log2ener: load energy from a log file
    - log2frc: load forces from a log file

    Args:
        project (str): the CP2K project name

    Returns:
        Dataset: a TIPS dataset

    """
    from os.path import exists
    from ase.data import atomic_numbers
    from tips.io.dataset import Dataset

    default_pattern = [
        ("pos", f"{project}-{cp2k_pos}", _index_xyz, _load_pos, _pos_spec),
        ("frc", f"{project}-{cp2k_frc}", _index_xyz, _load_frc, _frc_spec),
        ("cell", f"{project}-{cp2k_cell}", _index_cell, _load_cell, _cell_spec),
        ("ener", f"{project}-{cp2k_ener}", _index_ener, _load_ener, _ener_spec),
    ]

    indices, loaders, specs = {}, {}, {}
    for key, pattern, indexer, loader, spec in default_pattern:
        path = pattern.format(project=project)
        path = path if exists(path) else False
        if path:
            indices[key] = indexer(path)
            loaders[key] = loader
            specs.update(spec)

    def loader(i):
        data = {}
        for k in loaders.keys():
            data.update(loaders[k](indices[k][i]))
        return data


    sizes = {k:len(idx) for k,idx in indices.items()}
    try:
        assert len(set(sizes.values())) == 1, f"In consistent sizes {sizes}"
        files = list(sizes.keys())
        size = list(sizes.values())[0]
        logger.info(f"Indexed CP2K project {project}, size: {size}, files: {files}")
    except AssertionError as err:
        logger.exception(err)

    meta = {
        "fmt": "CP2K output",
        "size": size,
        "elem": set(loader(0)["elem"]),
        "spec": specs,
    }

    return Dataset(meta=meta, indexer=loader)
