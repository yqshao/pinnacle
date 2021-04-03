#!/usr/bin/env python3

from ase.io import iread
import contextlib
import io
import sys

def read(filename, format='auto', log=None, emap=None, units='real'):
    if filename.endswith('.yml'): #pinn handels yaml
        ds = pinn_reader(filename)
    else: #ase handels the rest
        kwargs = {}
        if filename.endswith('dump'): kwargs['units']=units
        ds = map(atoms2dict, iread(filename, **kwargs))
    if emap is not None:
        emap = {int(s.split(':')[0]):int(s.split(':')[1]) for s in emap.split(',')}
        map_elems = lambda data: dict(data, elems=[emap[e] for e in data['elems']])
        ds = map(map_elems, ds)
    if log is not None:
        update_energy = lambda struct, energy: dict(struct, e_data=energy)
        energies = read_lammps_log(log, units=units )
        ds = map(update_energy, ds, energies)
    return ds

def read_lammps_log(lammpsLog, units):
    """"Read lammps log file, for now we collect only TotEng"""
    from ase.calculators.lammps import convert
    with open(lammpsLog) as f:
        for line in f:
            if 'TotEng' in line:
                idx = line.split().index('TotEng')
                break
        for line in f:
            if not line.startswith('Loop'):
                yield convert(float(line.split()[idx]), 'energy', units, 'ASE')
            else:
                break

def pinn_reader(fname):
    """Load tfrecord dataset.

    Args:
       fname (str): filename of the .yml metadata file to be loaded.
       dtypes (dict): dtype of dataset.
    """
    import yaml, os
    import tensorflow as tf
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    with open(fname, 'r') as f:
        format_dict = (yaml.safe_load(f)['format'])
    dtypes = {k: format_dict[k]['dtype'] for k in format_dict.keys()}
    shapes = {k: format_dict[k]['shape'] for k in format_dict.keys()}
    feature_dict = {k: tf.io.FixedLenFeature([], tf.string) for k in dtypes}
    def parser(example): return tf.io.parse_single_example(example, feature_dict)
    def converter(tensors):
        tensors = {k: tf.io.parse_tensor(v, dtypes[k])
                   for k, v in tensors.items()}
        [v.set_shape(shapes[k]) for k, v in tensors.items()]
        return tensors
    tfr = '.'.join(fname.split('.')[:-1]+['tfr'])
    dataset = tf.data.TFRecordDataset(tfr).map(parser).map(converter)
    return dataset.as_numpy_iterator()


def get_writer(filename, format='pinn', **kwargs):
    if format == 'pinn':
        spec = {'elems':  {'shape':[None,],  'dtype':'int32'  },
                'coord':  {'shape':[None,3], 'dtype':'float32'},
                'cell':   {'shape':[3,3],    'dtype':'float32'},
                'e_data': {'shape':[],       'dtype':'float32'},
                'f_data': {'shape':[None,3], 'dtype':'float32'}}
        return pinn_writer(f'{filename}.yml', spec)
    elif format == 'lammps':
        return lammps_writer(f'{filename}.dump')
    elif format == 'xyz' or format == 'extxyz':
        return extxyz_writer(f'{filename}.xyz')
    else:
        raise NotImplementedError


class extxyz_writer():
    def __init__(self, fname):
        self.f = open(fname, 'w')

    def add(self, data):
        from ase import Atoms
        from ase.io import write
        from ase.calculators.singlepoint import SinglePointCalculator
        atoms = Atoms(data['elems'], cell=data['cell'], positions=data['coord'])
        atoms.calc = SinglePointCalculator(
            atoms=atoms,
            energy=data['e_data'],
            forces=data['f_data'])
        atoms.pbc=True
        write(self.f, atoms, format='extxyz', append='True')

    def finalize(self):
        self.f.close()


class lammps_writer():
    def __init__(self, fname):
        self.f = open(fname, 'w')
        self.count = 0

    def add(self, data):
        self.count+=1
        self.f.write(f"""ITEM: TIMESTEP
{self.count}
ITEM: NUMBER OF ATOMS
{len(data['elems'])}
ITEM: BOX BOUNDS pp pp pp
0.0000000000000000e+00 {data['cell'][0][0]:.16e}
0.0000000000000000e+00 {data['cell'][0][0]:.16e}
0.0000000000000000e+00 {data['cell'][0][0]:.16e}
ITEM: ATOMS id type x y z
""")
        self.f.writelines(
            f"{i+1} {elem} {coord[0]} {coord[1]} {coord[2]}\n"
            for i, (elem,coord)
            in enumerate(zip(data['elems'], data['coord'])))
    def finalize(self):
        self.f.close()


class pinn_writer():
    def __init__(self, fname, spec):
        import os
        import tensorflow as tf
        self.fname = fname
        self.spec = spec
        self.count = 0
        tfr = '.'.join(fname.split('.')[:-1]+['tfr'])
        self.writer = tf.io.TFRecordWriter(tfr)
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

    def _bytes_feature(self, value):
        import tensorflow as tf
        """Returns a bytes_list from a string / byte."""
        return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))

    def _serialize(self, tensors):
        import tensorflow as tf
        return {k: tf.io.serialize_tensor(tf.constant(v,dtype=self.spec[k]['dtype']))
                   for k, v in tensors.items()}

    def add(self, data):
        import tensorflow as tf
        tensors = self._serialize(data)
        example = tf.train.Example(
            features=tf.train.Features(
                feature={key: self._bytes_feature(val.numpy())
                         for key, val in tensors.items()}))
        self.writer.write(example.SerializeToString())
        self.count += 1

    def finalize(self):
        import yaml
        info_dict = {'n_sample': self.count}
        with open(self.fname, 'w') as f:
            yaml.safe_dump({'format': self.spec, 'info': info_dict}, f)


def atoms2dict(atoms):
    atomsDict = {
        'elems': atoms.numbers,
        'coord': atoms.positions,
        'cell': atoms.cell[:],
        'f_data': atoms.get_forces(),
        'e_data': atoms.get_potential_energy()
    }
    return atomsDict
