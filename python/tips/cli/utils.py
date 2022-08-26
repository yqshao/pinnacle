# -*- coding: utf-8 -*-

import click
from .common import load_opts, subsample_opts, load_ds_with_opts

def _gen_cp2k_inp(fin, fout, datum):
    """helper function to insert a structure to CP2K input"""
    from ase.data import chemical_symbols as symbol
    coord = [f'  {symbol[e]} {x} {y} {z}' for e, (x,y,z) in zip(datum['elem'], datum['coord'])]
    cell = [f'  {v} {x} {y} {z}' for v, (x, y, z) in zip('ABC', datum['cell'])]
    subsys = ['&COORD']+coord+['&END COORD']+['&CELL']+cell+['&END CELL']
    lines = open(fin).readlines()
    for idx, line in enumerate(lines):
        if '&END SUBSYS' in line:
            indent = len(line) - len(line.lstrip())
            break
    subsys = [' '*(indent+2) + l + '\n' for l in subsys]
    lines = lines[:idx] + subsys + lines[idx:]
    with open(fout, 'w') as f:
        f.writelines(lines)


@click.command(name="mkcp2kinp", short_help="make CP2K inputs from datasets")
@click.argument("inp", nargs=1)
@click.argument("dataset", nargs=1)
@load_opts
@click.option("--idx", default=-1, help="index of sample")
@click.option("--subsample", is_flag=True, default=False, help="activate subsample mode")
@subsample_opts
def mkcp2kinp(
        inp, dataset,
        fmt, emap, # load_opts
        idx, subsample, # options for this command
        strategy, nsample, psample, sort_key # subsample_opts
):
    ds = load_ds_with_opts(dataset, fmt, emap)
    if not subsample:
        _gen_cp2k_inp(inp, 'cp2k.inp', ds[idx])
    else:
        idx, subds = ds.subsample(strategy, nsample, psample, sort_key)
        for i, datum in zip(idx, subds):
            _gen_cp2k_inp(inp, f'{i}.inp', datum)

