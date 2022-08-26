# -*- coding: utf-8 -*-

import click
from .common import load_opts, shuffle_opts, write_opts, subsample_opts, load_ds_with_opts
from .utils import mkcp2kinp


@click.command(name="convert", short_help="convert datasets")
@click.argument("dataset", nargs=-1)
@load_opts
@write_opts
@shuffle_opts
def convert(
        dataset, fmt, emap, # load_opts
        output, ofmt, # write_opts,
        shuffle, seed #  shuffle_opts
):
    if not dataset: return
    ds = load_ds_with_opts(dataset, fmt, emap)
    if shuffle:
        ds.shuffle(seed=seed)
    ds.convert(output, fmt=ofmt)


@click.command(name="subsample", short_help="subsample datasets")
@click.argument("dataset", nargs=1)
@load_opts
@write_opts
@subsample_opts
def subsample(
        dataset, fmt, emap, # load_opts
        output, ofmt, # write_opts,
        strategy, nsample, psample, sort_key # subsample_opts
):
    ds = load_ds_with_opts(dataset, fmt, emap)
    idx, subds= ds.subsample(strategy, nsample, psample, sort_key)
    subds.convert(output, fmt=ofmt)


@click.command()
def version():
    import tips
    click.echo(f"TIPS version: {tips.__version__}")
