#!/usr/bin/env python3
import tips
import click
from tips.io import read, get_writer

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group()
def main():
    """TIPS CLI - A command line tool for manipulating of AML data"""
    pass

@click.command()
def version():
    click.echo(f'TIPS version: {tips.__version__}')

@click.command(name='convert', context_settings=CONTEXT_SETTINGS, short_help='convert datasets')
@click.argument('filename', nargs=-1)
@click.option('-f', '--format', metavar='', default='auto', help='input format')
@click.option('-o', '--output', metavar='', default='dataset')
@click.option('-of', '--oformat', metavar='', default='pinn', help='output format')
@click.option('--shuffle', metavar='', default=True, help='shuffle the dataset')
@click.option('--seed', metavar='', default='0',  type=int, help='seed for random number (int)')
@click.option('--log', metavar='', default=None, help='lammps log (for energies)')
@click.option('--units', metavar='', default=None, help='lammps units')
@click.option('--emap', metavar='', default=None, help='remap lammps elements, e.g. "1:1,2:8"')
def convertds(filename, format, output, oformat, shuffle, seed, log, units, emap):
    import itertools, random, math
    if ':' in output:
        split = True
        writers = [get_writer(s.split(':')[0], format=oformat) for s in output.split(',')]
        weights = [float(s.split(':')[1]) for s in output.split(',')]
    else:
        split = False
        writers = [get_writer(output, format=oformat)]
        writerList = itertools.repeat(writers[0])

    for fname in filename:
        dataset = read(fname, format=format, log=log, emap=emap, units=units)
        if split:
            dataset, ds4count = itertools.tee(dataset)
            count = sum(1 for _ in ds4count)
            writerList = sum([[writer]*math.ceil(count*weight/sum(weights))
                              for writer, weight in zip(writers, weights)],[])
            if shuffle:
                random.seed(seed)
                random.shuffle(writerList)
        for datum, writer in zip(dataset, writerList):
            writer.add(datum)
    [writer.finalize() for writer in writers]


@click.command(name='filter', context_settings=CONTEXT_SETTINGS, short_help='filter datasets')
@click.argument('filename', nargs=-1)
@click.option('--log', metavar='', default=None, help='lammps log (for energies)')
@click.option('--units', metavar='', default=None, help='lammps units')
@click.option('--emap', metavar='', default=None, help='remap lammps elements, e.g. "1:1,2:8"')
@click.option('-f', '--format', metavar='', default='auto', help='input format')
@click.option('-o', '--output', metavar='', default='dataset')
@click.option('-of', '--oformat', metavar='', default='pinn', help='output format')
@click.option('-a', '--algo', metavar='', default='naive', help='filtering algorithm')
@click.option('-vmin', '--val-min', metavar='', default='', help='minimum value')
@click.option('-vmax', '--val-max', metavar='', default='', help='maximum value')
@click.option('-amin', '--abs-min', metavar='', default='', help='minimum absolute value')
@click.option('-amax', '--abs-max', metavar='', default='', help='maximum absolute value')
def filterds(filename, log, units, emap, format, output, oformat, algo, **kwargs):
    """\b
    Algorithms available:
    - 'naive': filter by the error tolerance
    - 'qbc': filter by the standard deviation across datasets
    """
    import itertools
    import numpy as np
    from tips.filters import qbc_filter, naive_filter

    writer = get_writer(output, format=oformat)
    datasets = [read(fname, log=log, format=format, emap=emap, units=units)
                for fname in filename]
    if algo=='qbc':
        ds = qbc_filter(zip(*datasets), **kwargs)
    elif algo=='naive':
        ds = naive_filter(itertools.chain(*datasets), **kwargs)
    else:
        raise f"Unknown filter {algo}"

    indices = []
    for idx, data in ds:
        indices.append(idx)
        writer.add(data)

    np.savetxt(f'{output}.idx', indices, fmt='%i')
    writer.finalize()

main.add_command(convertds)
main.add_command(filterds)
main.add_command(version)

if __name__ == '__main__':
    main()
