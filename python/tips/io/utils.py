# -*- coding: utf-8 -*-

import tips

def tips_loader(name):
    def register_loader(func):
        tips.io.loaders[name] = func
        return func
    return register_loader

def tips_convertor(name):
    def register_convertor(func):
        tips.io.convertors[name] = func
        return func
    return register_convertor


def list_loader(func):
    """Wraps a function such that a joined dataset is returned when a list
    inputs is provided"""
    from functools import wraps, reduce
    from collections.abc import Iterable

    @wraps(func)
    def wrapper(path, *args, **kwargs):
        if isinstance(path, Iterable):
            all_ds = [func(p, *args, **kwargs) for p in path]
            return reduce(lambda x, y: x.join(y), all_ds)
        else:
            return func(path)

    return wrapper
