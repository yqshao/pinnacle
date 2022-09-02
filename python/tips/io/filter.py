# -*- coding: utf-8 -*-

import re
import numpy as np

def filters2fn(filters):
    """
    Compiles a lisf of filter flags into a function for datum

    Filters should follow one of the following pattterns, without space:

    - [key][compare][tol], i.e. energy>-10
    - [op]([key])[compare][tol], i.e. abs(force)>0.1

    All filters are joined with `all()` at this point, and when the target
    array is not a scalar, the `np.all` will be applied.
    """
    _comparers = {
        '>': lambda x, y: np.all(x>y),
        '=': lambda x, y: np.all(x==y),
        '<': lambda x, y: np.all(x<y),
    }
    _ops = {
        'abs': lambda x, k: np.abs(x[k]),
        'peratom': lambda x, k: x[k]/np.shape(x['elem'])[0],
    }
    _filters = []

    regex = r"^(?:([a-zA-Z-_]+)|(abs|peratom)\(([a-zA-Z-_]+)\))([>=<])([+-]?[0-9]*[.]?[0-9]+)$"
    # matches [          1:tag or       2:op,            3:tag][4:compare]         [5:float]
    for filter in filters:
        match = re.match(regex, filter)
        assert match, f"filter {filter} not recognized"
        if match.group(1):
            op = lambda x: x
            tag = match.group(1)
            comp = _comparers[match.group(4)]
            tol = float(match.group(5))
        else:
            op = _ops[match.group(2)]
            tag = match.group(3)
            comp = _comparers[match.group(4)]
            tol = float(match.group(5))
        _filter = lambda x, comp=comp, op=op, tol=tol, tag=tag:\
            comp(op(x, tag), tol) # ^ making  _filter proper closures
        _filters.append(_filter)
    filter_fn = lambda datum: all(_filter(datum) for _filter in _filters)
    return filter_fn
