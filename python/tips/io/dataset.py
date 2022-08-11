# -*- coding: utf-8 -*-

"""This module implements the Dataset class, for data manuputaion in TIPS

"""
import tips
from copy import deepcopy


class Dataset:
    def __init__(self, generator=None, indexer=None, meta=None):
        """Dataset should be initialized with a generator, optionally, an
        inedexer

        """
        self.meta = meta
        self.indexer = indexer
        if indexer is not None and generator is None:

            def generator():
                for i in range(self.meta["size"]):
                    yield self.indexer(i)

        self.generator = generator

    def __getitem__(self, key):
        if self.indexer is None:
            raise NotImplementedError("The dataset is not indexible")

        if isinstance(key, slice):
            keys = list(range(self.meta["size"]))[key]
            meta = deepcopy(self.meta)
            meta["size"] = len(keys)
            meta["fmt"] = "TIPS sliced"

            def new_indexer(i):
                return self.indexer(keys[i])

            return Dataset(indexer=new_indexer, meta=meta)
        elif isinstance(key, int):
            if key < 0:
                key = key + self.meta["size"]
            if key >= self.meta["size"] or key < 0:
                raise IndexError("Index out of range.")
            return self.indexer(key)
        else:
            raise NotImplementedError(f"Cannot index a dataset with {type(key)}")

    def __iter__(self):
        return self.generator()

    def __repr__(self):
        from ase.data import chemical_symbols

        spec_repr = "\n     ".join(
            f"{k}: {v['shape']}, {v['dtype']}" for k, v in self.meta["spec"].items()
        )
        repr = f"""<tips.io.Dataset
 fmt: {self.meta['fmt']}
 size: {self.meta['size']}
 elem: {', '.join(chemical_symbols[i] for i in sorted(self.meta['elem']))}
 spec:
     {spec_repr}>"""
        return repr

    def convert(self, *args, fmt=None, **kwargs):
        """Convert the dataset to a known format the format, the format writer
        should be one of the registered writers, with the `tips.convertor`
        decorator. See `docs/python/io/#custom-readerwriter` for an example.
        Args:
           fmt: target format as defined in tips.io.convertors

        Returns:
           dependent on the convertor
        """
        convertor = tips.io.convertors[fmt]
        return convertor(self, *args, **kwargs)

    def split(self, splits, shuffle=True, seed=0):
        """Split the dataset into randomly splitted formats

        Args:
            splits (dict): dictionary of fractions
            shuffle (bool): shuffle the dataset
            seed (int): random seed for splitting

        Return:
            dict of Datasets
        """
        return splitted

    def join(self, ds):
        """Joins two datasets

        Args:
            ds: a tips dataset object

        Return:
            joined dataset
        """

        # check the compatibility of two datsets
        assert set(self.meta["spec"].keys()) == set(self.meta["spec"].keys())
        meta = deepcopy(self.meta)
        meta["fmt"] = "TIPS joined"
        meta["size"] = self.meta["size"] + ds.meta["size"]
        meta["elem"] = self.meta["elem"].union(ds.meta["elem"])
        if self.indexer is not None and ds.indexer is not None:

            def indexer(i):
                if i < self.meta["size"]:
                    return self.indexer(i)
                else:
                    return ds.indexer(i - self.meta["size"])

            return Dataset(indexer=indexer, meta=meta)

    def filter(self, func):
        """Filter the dataset with a boolean function

        Args:
            func: a function that takes a datum and returns a bool

        Return:
            filtered dataset
        """
        assert set(self.meta["spec"].keys()) == set(self.meta["spec"].keys())
        meta = deepcopy(self.meta)
        meta["fmt"] = "TIPS filtered"
        meta["size"] = None

        def generator():
            for datum in self.generator():
                if func(datum):
                    yield datum

        return Dataset(generator=generator, meta=meta)

    def map_elems(self, emap):
        """Maps the elements of a dataset according to some typing rule

        Args:
            emap: a dict or a LAMMPS data file

        Return:
            mapped dataset
        """
        import numpy as np
        from ase.data import atomic_masses

        if isinstance(emap, dict):
            edict = emap

        edict = {}
        with open(emap) as f:
            while True:
                l = f.readline()
                if l.startswith("Masses"):
                    break
            f.readline()
            while True:
                line = f.readline().split()
                if len(line) == 0:
                    break
                else:
                    idx = int(line[0])
                    mass = float(line[1])
                    edict[idx] = max(1, np.argmin(np.abs(atomic_masses - mass)))
        def mapper(datum):
            datum = datum.copy()
            datum["elem"] = [edict[i] for i in datum["elem"]]
            return datum
        meta = deepcopy(self.meta)
        meta["fmt"] = "TIPS mapped"
        meta["elem"] = set([edict[i] for i in meta["elem"]])
        if self.indexer is None:
            def generator():
                for datum in self.generator():
                    yield mapper(datum)
            return Dataset(generator=generator, meta=meta)
        else:
            def indexer(i):
                return mapper(self.indexer(i))
            return Dataset(indexer=indexer, meta=meta)
