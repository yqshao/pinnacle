#!/usr/bin/env python3

def read(filename, format='auto', log=None, emap=None, top=None, units='real'):
    from ase.io import iread
    if format=='auto':
        if filename.endswith('.yml'):
            format='pinn'
        elif filename.endswith('.dump'):
            format='dump'
        elif filename.endswith('.trr'):
            format='trr'

    if format=='pinn':
        ds = pinn_reader(filename)
    elif format=='dump':
        ds = map(atoms2dict, iread(filename, units=units))
        if log is not None:
            update_energy = lambda struct, energy: dict(struct, e_data=energy)
            energies = read_lammps_log(log, units=units )
            ds = map(update_energy, ds, energies)
    elif format=='trr':
        ds = read_trr(filename, top)
        if log is not None:
            update_energy = lambda struct, energy: dict(struct, e_data=energy)
            energies = read_gromacs_xvg(log)
            ds = map(update_energy, ds, energies)
    else:
        ds = map(atoms2dict, iread(filename))

    if emap is not None:
        emap = {int(s.split(':')[0]):int(s.split(':')[1]) for s in emap.split(',')}
        map_elems = lambda data: dict(data, elems=[emap[e] for e in data['elems']])
        ds = map(map_elems, ds)
    return ds


def get_writer(filename, format='xyz', **kwargs):
    if format == 'pinn':
        spec = {'elems':  {'shape':[None,],  'dtype':'int32'  },
                'coord':  {'shape':[None,3], 'dtype':'float32'},
                'cell':   {'shape':[3,3],    'dtype':'float32'},
                'e_data': {'shape':[],       'dtype':'float32'},
                'f_data': {'shape':[None,3], 'dtype':'float32'}}
        return pinn_writer(f'{filename}.yml', spec)
    elif format == 'dump':
        return lammps_writer(f'{filename}.dump')
    else:
        return ase_writer(f'{filename}.{format}')


def read_trr(filename, top):
    import MDAnalysis as mda
    from ase.data import atomic_numbers
    u = mda.Universe(top, filename)
    elems = [atomic_numbers[e] for e in u.atoms.types]
    for ts in u.trajectory:
        data = {
            'cell':   ts.triclinic_dimensions,
            'elems':  elems,
            'coord':  ts.positions,
            'f_data': ts.forces,
            'e_data': 0.0,
        }
        yield data


def read_gromacs_xvg(fname):
    import numpy as np
    energies = np.loadtxt(fname, comments=['#','@'])[:,1]
    for energy in energies:
        yield energy # should be in kcal/mol


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


class ase_writer():
    def __init__(self, fname):
        self.fname = fname
        self.traj = []

    def add(self, data):
        from ase import Atoms
        from ase.calculators.singlepoint import SinglePointCalculator
        atoms = Atoms(data['elems'], cell=data['cell'], positions=data['coord'])
        atoms.calc = SinglePointCalculator(
            atoms=atoms,
            energy=data['e_data'],
            forces=data['f_data'])
        atoms.pbc=True
        self.traj.append(atoms)

    def finalize(self):
        from ase.io import write
        write(self.fname, self.traj)


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
    import numpy as np
    try:
        energy = atoms.get_potential_energy()
    except:
        energy = 0.0
    try:
        forces = atoms.get_forces()
    except:
        forces = np.zeros_like(atoms.positions)

    atomsDict = {
        'elems': atoms.numbers,
        'coord': atoms.positions,
        'cell': atoms.cell[:],
        'f_data': forces,
        'e_data': energy,
    }
    return atomsDict
