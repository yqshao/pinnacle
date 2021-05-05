#!/usr/bin/env python3
import numpy as np

filter_fns = {
    'val_max': lambda data, tol: np.any(data>tol),
    'val_min': lambda data, tol: np.any(data<tol),
    'abs_max': lambda data, tol: np.any(np.abs(data)>tol),
    'abs_min': lambda data, tol: np.any(np.abs(data)<tol)}


def get_filters(data_parser, **kwargs):
    """ Build a list of filter functions"""
    filters = []
    for key, val in kwargs.items():
        if not val:
            continue
        for tag, tol in map(lambda x: x.split(':'), val.split(',')):
            if tag.startswith('!'):
                filter_fn = lambda data, tol: not filter_fns[key](data, tol)
                tag = tag[1:]
            else:
                filter_fn = filter_fns[key]
            fn = lambda data, tag=tag, tol=tol, parser=data_parser, filter_fn=filter_fn:\
                filter_fn(parser(data, f'{tag}_data'), float(tol))
            filters.append(fn)
    return filters


def qbc_filter(ds, **kwargs):
    data_parser = lambda data, tag: np.std([d[tag] for d in data], axis=0)
    filters = get_filters(data_parser, **kwargs)
    for i, data in enumerate(ds):
        assert len(data)>=2
        assert np.all([np.allclose(data[0]['coord'], d['coord']) for d in data[1:]])
        if not np.any([fn(data) for fn in filters]):
            yield i, data[0]


def naive_filter(ds, **kwargs):
    data_parser = lambda data, key: data[key]
    filters = get_filters(data_parser, **kwargs)
    for i, data in enumerate(ds):
        print([fn(data) for fn in filters])
        if not np.any([fn(data) for fn in filters]):
            yield i, data
