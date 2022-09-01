# -*- coding: utf-8 -*-
#
"""This module implements the loader for CP2K log files"""

import numpy as np
from tips.io.utils import list_loader
from ase.units import create_units

# This is to follow the CP2K standard to use CODATA 2006, which differs from the
# the defaults of ASE (as of ASE ver 3.23 and CP2K v2022.1, Sep 2022)
units = create_units('2006')

def _index_pattern(fname, pattern):
    import mmap, re
    f = open(fname, 'r')
    m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
    locs = [match.span()[0] for match in
            re.finditer(pattern, m)]
    return locs


def _index_energy(fname):
    import mmap, re
    f = open(fname, 'r')
    regex = r'ENERGY\|\ Total FORCE_EVAL.*:\s*([-+]?\d*\.?\d*)'
    energies = [float(e)*units['Hartree'] for e in re.findall(regex, f.read())]
    f.close()
    return energies

def _load_cell(fname, loc):
    f = open(fname, 'r')
    cell = []
    f.seek(loc)
    assert (
        f.readline().startswith(' CELL| Volume [angstrom^3]:')
    ), "Unknown format of CP2K log, aborting"
    for i in range(3):
        l = f.readline().strip()
        cell.append(l.split()[4:7])
    return {'cell': np.array(cell, np.float)}

def _load_coord(fname, loc):
    f = open(fname, 'r')
    coord, elem = [], []
    f.seek(loc)
    assert (
        f.readline().startswith(' MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom') &
        f.readline().startswith('\n')  &
        f.readline().startswith(' Atom  Kind  Element')
    ), "Unknown format of CP2K log, aborting"
    l = f.readline().strip()
    while l:
        coord.append(l.split()[4:7])
        elem.append(l.split()[3])
        l = f.readline().strip()
    f.close()
    return {'coord': np.array(coord, np.float), 'elem': np.array(elem, np.int)}


def _load_force(fname, loc):
    f = open(fname, 'r')
    data = []
    f.seek(loc)
    assert (
        f.readline().startswith(' ATOMIC FORCES in [a.u.]') &
        f.readline().startswith('\n')  &
        f.readline().startswith(' # Atom   Kind   Element')
    ), "Unknown format of CP2K log, aborting"
    l = f.readline().strip()
    while not l.startswith('SUM OF'):
        data.append(l.split()[3:])
        l = f.readline().strip()
    f.close()
    return {'force': np.array(data, np.float)*units['Hartree']/units['Bohr']}


def _load_stress(fname, loc):
    f = open(fname, 'r')
    data = []
    f.seek(loc)
    assert (
        f.readline().startswith(' STRESS| Analytical stress tensor [GPa]') &
        f.readline().startswith(' STRESS|                        x')
    ), "Unknown format of CP2K log, aborting"
    for i in range(3):
        l = f.readline().strip()
        data.append(l.split()[2:])
    f.close()
    return {'stress': -np.array(data, np.float)*units['GPa']}

@list_loader
def load_cp2klog(fname):
    """Loads cp2k-formatted logs

    Args:
        fname (str): the path to CP2K log file

    Returns:
        Dataset: a TIPS dataset

    """
    from os.path import exists
    from ase.data import atomic_numbers
    from tips.io.dataset import Dataset

    specs = {
        'cell': {"cell": {"dtype": "float", "shape": [3, 3]}},
        'coord': {
            "elem": {"dtype": "int", "shape": [None]},
            "coord": {"dtype": "float", "shape": [None, 3]},
        },
        'energy': {"energy": {"dtype": "float", "shape": []}},
        'force': {"force": {"dtype": "float", "shape": [None, 3]}},
        'stress': {"stress": {"dtype": "float", "shape": [3, 3]}}
    }

    indexers = {
        'cell': lambda fname: _index_pattern(fname, b' CELL\| Volume \[angstrom\^3\]'),
        'coord':lambda fname: _index_pattern(fname, b' MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom'),
        'energy': _index_energy,
        'force': lambda fname: _index_pattern(fname, b' ATOMIC FORCES in \[a.u.\]'),
        'stress': lambda fname: _index_pattern(fname, b' STRESS\| Analytical stress tensor \[GPa\]'),
    }

    indices, spec = {}, {}
    for k, v in indexers.items():
        idx =  v(fname)
        if len(idx) == 0:
            continue
        else:
            indices[k] = idx
            spec.update(specs[k])

    sizes = {k:len(idx) for k,idx in indices.items()}
    assert len(set(sizes.values())) == 1, f"Inconsistent sizes {sizes}"
    size = list(sizes.values())[0]
    keys = list(sizes.keys())

    loaders = {
        'cell': _load_cell,
        'coord': _load_coord,
        'energy': lambda fname, energy: {'energy': energy},
        'force': _load_force,
        'stress': _load_stress,
    }

    def loader(i):
        data = {}
        for k in loaders.keys():
            data.update(loaders[k](fname, indices[k][i]))
        return data

    meta = {
        "fmt": "CP2K log",
        "size": size,
        "spec": spec,
    }

    return Dataset(meta=meta, indexer=loader)
