# -*- coding: utf-8 -*-


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
    for k in spec.keys():
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
    info_dict = {
        "n_sample": ds.meta["size"],
        "elems": [int(e) for e in list(ds.meta["elem"])],
    }
    with open(f"{fname}.yml", "w") as f:
        yaml.safe_dump({"format": spec, "info": info_dict}, f)
