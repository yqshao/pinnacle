# -*- coding: utf-8 -*-
from .utils import list_loader

_pinn2tips = {
    'cell': 'cell',
    'coord': 'coord',
    'elems': 'elem',
    'e_data': 'energy',
    'f_data': 'force',
    's_data': 'stress',
}


@list_loader
def load_pinn_tfr(fname):
    import sys, yaml
    import numpy as np
    import tensorflow as tf
    from tensorflow.python.lib.io.file_io import FileIO
    from tips.io.dataset import Dataset
    with FileIO(fname, 'r') as f:
        ds_spec = yaml.safe_load(f)
    size = ds_spec['info']['n_sample']
    format_dict = ds_spec['format']
    dtypes = {k: format_dict[k]['dtype'] for k in format_dict.keys()}
    shapes = {k: format_dict[k]['shape'] for k in format_dict.keys()}
    feature_dict = {k: tf.io.FixedLenFeature([], tf.string) for k in dtypes}
    spec = {_pinn2tips[k]:v for k,v in format_dict.items()}
    meta = {'fmt':'PiNN-style TFRecord', 'size':size, 'spec':spec}
    def parser(example):
        return tf.io.parse_single_example(example, feature_dict)
    def converter(tensors):
        tensors = {k: tf.io.parse_tensor(v, dtypes[k])
                   for k, v in tensors.items()}
        [v.set_shape(shapes[k]) for k, v in tensors.items()]
        tensors = {_pinn2tips[k]:v for k,v in tensors.items()}
        return tensors
    tfr = '.'.join(fname.split('.')[:-1]+['tfr'])
    dataset = tf.data.TFRecordDataset(tfr).map(parser).map(converter)
    generator = dataset.as_numpy_iterator
    return Dataset(meta=meta, generator=generator)


def _bytes_feature(value):
    """Returns a bytes_list from a string / byte."""
    import tensorflow as tf

    return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))


def _serialize(tensors, spec):
    import tensorflow as tf

    return {
        k: tf.io.serialize_tensor(tf.constant(v, dtype=spec[k]["dtype"]))
        for k, v in tensors.items()
    }


def ds2pinn(ds, fname):
    """Converts a dataset"""
    import yaml
    import tensorflow as tf
    from copy import deepcopy

    tfr = f"{fname}.tfr"
    writer = tf.io.TFRecordWriter(tfr)

    def tips2pinn(dict_tips):
        convert = {
            "force": "f_data",
            "energy": "e_data",
            "stress": "s_data",
            "elem": "elems",
        }
        dict_pinn = {
            k if (k not in convert) else convert[k]: v for k, v in dict_tips.items()
        }
        return dict_pinn

    spec = deepcopy(tips2pinn(ds.meta["spec"]))
    for k in spec.keys(): # default to single precision numbers
        if '32' not in spec[k]["dtype"] and '64' not in spec[k]["dtype"]:
            spec[k]["dtype"] += "32"

    for idx, data in enumerate(ds):
        tensors = _serialize(tips2pinn(data), spec)
        example = tf.train.Example(
            features=tf.train.Features(
                feature={
                    key: _bytes_feature(val.numpy()) for key, val in tensors.items()
                }
            )
        )
        writer.write(example.SerializeToString())

    # finalize
    info_dict = {"n_sample": idx+1}
    if 'elem' in ds.meta:
        info_dict["elems"]: [int(e) for e in list(ds.meta["elem"])]
    with open(f"{fname}.yml", "w") as f:
        yaml.safe_dump({"format": spec, "info": info_dict}, f)
