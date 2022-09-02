#!/usr/bin/env python

import click

def _group_options(*options):
    def wrapper(function):
        for option in reversed(options):
            function = option(function)
        return function
    return wrapper

load_opts = _group_options(
    click.option("-f", "--fmt", default="auto", help="input format"),
    click.option("-em", "--emap", default=None, help="lammps data for elem"),
)

def load_ds_with_opts(dataset, fmt, emap):
    from tips.io import load_ds
    ds = load_ds(dataset, fmt=fmt)
    if emap is not None:
        ds = ds.map_elems(emap)
    return ds

shuffle_opts = _group_options(
    click.option("--shuffle", is_flag=True, default=False, help="shuffle the dataset"),
    click.option("--seed", default=0, help="seed for shuffling"),
)

write_opts = _group_options(
    click.option("-o", "--output" , default="output", help="output file"),
    click.option("-of", "--ofmt", default="extxyz", help="output format"),
)

subsample_opts = _group_options(
    click.option('--strategy', default='uniform', help="subsampling strategy"),
    click.option('--nsample', default=None, type=int, help="number to subsample"),
    click.option('--psample', default=None, type=float, help="percentage to subsample"),
    click.option('--sort-key', default='force_std', help="key for the sorted strategy"),
)

filter_opts = _group_options(
    click.option('--filter', "filters", default=[], multiple=True, help="dataset filters"),
)
