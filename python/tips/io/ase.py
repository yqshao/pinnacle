# -*- coding: utf-8 -*-

"""Data Loader for ASE Trajectory Objects and .traj files

The readers and writers are patched from the original ASE implementation such
that extra atomic and structrural featurse can be written.

"""

import logging
import numpy as np
from mock import patch
from tips.io.utils import list_loader
from ase.calculators.calculator import all_properties

logger = logging.getLogger('tips')

# patches for ASE IO modules
extra_properties = [
    f'{prop}_{extra}'
    for prop in ['energy', 'forces', 'stress']
    for extra in ['avg','std','bias']
]
per_atom_properties = [f'{prop}{suffix}'
                       for prop in ['forces']
                       for suffix in ['', '_avg','_std','_bias']]
per_config_properties = [f'{prop}{suffix}'
                         for prop in ['energy', 'stress']
                         for suffix in ['', '_avg','_std','_bias']]
all_prop_patch = patch('ase.calculators.calculator.all_properties', all_properties+extra_properties)
atom_prop_patch = patch('ase.io.extxyz.per_atom_properties', per_atom_properties)
struc_prop_patch = patch('ase.io.extxyz.per_config_properties', per_config_properties)


def _tips2ase(k):
    """naming translator form tips to ase"""
    if 'force' in k:
        return k.replace('force', 'forces')
    else:
        return k

def _ase2tips(k):
    """naming translator form tips to ase"""
    if 'forces' in k:
        return k.replace('forces', 'force')
    else:
        return k


def load_ase(traj, skim=True, atomic=['force']):
    """Loading Dataset from a ASE trajectory

    Args:
        traj: a trajectory object (list of atoms)
        skim: skim the trajectory to get the element information
        atomic: list of keys to match the atomic features

    Return:
        A TIPS Dataset
    """
    from tips.io.dataset import Dataset
    traj = traj.copy()

    spec = {
        "elem": {"shape": [None], "dtype": "int"},
        "coord": {"shape": [None, 3], "dtype": "float"},
    }

    if traj[0].pbc.all(): # only adds cell for periodic systemcs
        spec['cell'] =  {"shape": [3, 3], "dtype": "float"}

    k_atomic, k_structural = [], []
    for k, v in traj[0].calc.results.items():
        # first check the dtype with numpy
        v=np.array(v)
        if np.issubdtype(v.dtype, np.integer):
            dtype = 'int'
        elif np.issubdtype(v.dtype, np.float):
            dtype = 'float'
        else:
            logger.info(f'Skipping unrecognizable result {k}')
            break
        # then determine the shape
        if any([(ka in k) for ka in atomic]):
            shape = [None, *np.shape(v)[1:]]
            k_atomic.append(k)
        else:
            shape = list(np.shape(v))
            k_structural.append(k)
        spec[_ase2tips(k)] = {'shape':shape, 'dtype':dtype}
    if k_atomic:
        logger.info(f'Inferring atomic features: {k_atomic}')
    if k_structural:
        logger.info(f'Inferring structural features: {k_structural}')

    meta = {
        "fmt": "ASE Dataset",
        "size": len(traj),
        "spec": spec
    }

    def indexer(i):
        atoms = traj[i]
        datum = {
            "elem": atoms.numbers,
            "coord": atoms.positions,
        }
        if "cell" in meta["spec"]:
            datum["cell"] = atoms.cell[:]
        for k, v in meta["spec"].items():
            if k in ["cell", "elem", "coord"]:
                continue # hard-coded features
            datum[k] = atoms.calc.results[_tips2ase(k)]
        return datum

    ds = Dataset(meta=meta, indexer=indexer)
    if skim:
        ds.skim()
    return ds


@list_loader
def load_asetraj(fname, index=':', **kwargs):
    """Loading Dataset from a ASE trajectory

    Args:
        fname: a trajectory file
        index: indices of the trajectory to load
        **kwargs: arguments applicable to load_ase

    Return:
        a TIPS Dataset
    """
    from ase import Atoms
    with all_prop_patch:
        from ase.io import read
        traj = read(fname, index=index)
        if isinstance(traj, Atoms):
            traj = [traj]
    return load_ase(traj, **kwargs)


def ds2ase(dataset):
    from ase import Atoms
    from ase.calculators.singlepoint import SinglePointCalculator
    traj = []
    for datum in dataset:
        atoms = Atoms(datum['elem'], positions=datum['coord'])
        if 'cell' in datum:
            atoms.cell = datum['cell']
            atoms.pbc = True
        results = {_tips2ase(k):v for k,v in datum.items()
                   if k not in ['elem', 'coord', 'cell']}
        calc = SinglePointCalculator(atoms, **results)
        atoms.calc = calc
        traj.append(atoms)
    return traj


def ds2asetraj(dataset, fname):
    """Writes a dataset to an ASE trajectory file, with extra columns

    Args:
        dataset: a TIPS Datset object
        fname: output file name ('.traj' will be appended if not alreay there)
    """
    from ase.calculators.calculator import all_properties
    if not fname.endswith('.traj'):
        fname += '.traj'
    extra_properties = [
        f'{prop}_{extra}'
        for prop in ['energy', 'forces', 'stress']
        for extra in ['avg','std','bias']
    ]
    traj = ds2ase(dataset)
    with all_prop_patch:
        from ase.io import write
        write(fname, traj)


def ds2extxyz(dataset, fname):
    """Writes a dataset to an extended-xyz file, with extra columns

    Args:
        dataset: a TIPS Datset object
        fname: output file name ('.xyz' will be appended if not alreay there)
    """
    if not fname.endswith('.xyz'):
        fname += '.xyz'
    traj = ds2ase(dataset)
    for atoms in traj:
        results = atoms.calc.results
        if 'stress' in results and results['stress'].shape == (3, 3):
            results['stress'] = results['stress'][[0,1,2,1,0,0],[0,1,2,2,2,1]]
    with all_prop_patch, atom_prop_patch, struc_prop_patch:
        from ase.io import write
        write(fname, traj)
