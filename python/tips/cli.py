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
@click.option('--log', metavar='', default=None, help='lammps log (for energies)')
@click.option('--units', metavar='', default=None, help='lammps units')
@click.option('--emap', metavar='', default=None, help='remap lammps elements, e.g. "1:1,2:8"')
@click.option('-f', '--format', metavar='', default='auto', help='input format')
@click.option('-o', '--output', metavar='', default='dataset')
@click.option('-of', '--oformat', metavar='', default='pinn', help='output format')
def convertds(filename, log, units, format, output, oformat, emap):
    writer = get_writer(output, format=oformat)
    for fname in filename:
        dataset = read(fname, format=format, log=log, emap=emap, units=units)
        with click.progressbar(dataset, show_pos=True, bar_template='Converting: %(info)s structures.') as ds:
            for datum in ds:
                writer.add(datum)
    writer.finalize()

@click.command(name='split', context_settings=CONTEXT_SETTINGS, short_help='split datasets')
@click.argument('filename', nargs=-1)
@click.option('--log', metavar='', default=None, help='lammps log (for energies)')
@click.option('--units', metavar='', default=None, help='lammps units')
@click.option('--emap', metavar='', default=None, help='remap lammps elements, e.g. "1:1,2:8"')
@click.option('-s', '--splits', metavar='', default='train:8,test:2', help='name and ratio of splits')
@click.option('--shuffle', metavar='', default=True, help='shuffle the dataset')
@click.option('--seed', metavar='', default='0',  type=int, help='seed for random number (int)')
@click.option('-f', '--format', metavar='', default='auto', help='input format')
@click.option('-of', '--oformat', metavar='', default='pinn', help='output format')
def splitds(filename, log, units, format, splits, shuffle, seed, oformat, emap):
    import random, itertools, math, time
    writers = [get_writer(s.split(':')[0], format=oformat) for s in splits.split(',')]
    weights = [float(s.split(':')[1]) for s in splits.split(',')]
    for fname in filename:
        dataset = read(fname, format=format, log=log, emap=emap, units=units)
        dataset, ds4count = itertools.tee(dataset)
        count = sum(1 for _ in ds4count)
        writerList = sum([[writer]*math.ceil(count*weight/sum(weights))
                          for writer, weight in zip(writers, weights)],[])
        if shuffle:
            random.seed(seed)
            random.shuffle(writerList)
        with click.progressbar(dataset, length=count, show_pos=True) as ds:
            for datum, writer in zip(ds, writerList):
                writer.add(datum)
    [writer.finalize() for writer in writers]

@click.command(name='qbc', context_settings=CONTEXT_SETTINGS, short_help='query by committee')
@click.argument('filename', nargs=-1)
@click.option('--log', metavar='', default=None, help='lammps log (for energies)')
@click.option('--units', metavar='', default=None, help='lammps units')
@click.option('--emap', metavar='', default=None, help='remap lammps elements, e.g. "1:1,2:8"')
@click.option('-o', '--output', metavar='', default='dataset')
@click.option('-f', '--format', metavar='', default='auto', help='input format')
@click.option('-of', '--oformat', metavar='', default='pinn', help='output format')
@click.option('-t', '--tags', metavar='', default='e,f', help='tags to compute variance')
def qbc(filename, log, units, emap, format, output, oformat, tags):
    """Label errors according to variable across datasets

    Not that energies for lammps dumps will not be handeled correctly for now
    """
    import numpy as np
    writer = get_writer(output, format=oformat)
    ds = [read(fname, format=format, emap=emap, units=units) for fname in filename]
    for data in zip(*ds):
        for tag in tags.split(','):
            for other in data[1:]:
                assert np.allclose(other['coord'],data[0]['coord']), 'Inputs does not match'
            error = np.std([datum[f'{tag}_data'] for datum in data], axis=0)
            data[0].update({f'{tag}_data': error})
        writer.add(data[0])
    writer.finalize()


@click.command(name='filter', context_settings=CONTEXT_SETTINGS, short_help='filter datasets')
@click.argument('filename', nargs=-1)
@click.option('--log', metavar='', default=None, help='lammps log (for energies)')
@click.option('--units', metavar='', default=None, help='lammps units')
@click.option('--emap', metavar='', default=None, help='remap lammps elements, e.g. "1:1,2:8"')
@click.option('-f', '--format', metavar='', default='auto', help='input format')
@click.option('-o', '--output', metavar='', default='dataset')
@click.option('-of', '--oformat', metavar='', default='pinn', help='output format')
@click.option('-vmin', '--val-min', metavar='', default='', help='minimum value')
@click.option('-vmax', '--val-max', metavar='', default='', help='maximum value')
@click.option('-amin', '--abs-min', metavar='', default='', help='minimum absolute value')
@click.option('-amax', '--abs-max', metavar='', default='', help='maximum absolute value')
def filterds(filename, log, units, emap, format, output, oformat, **kwargs):
    """\b
    Algorithms available:
    - 'naive': filter by the error tolerance
    - 'fps': furthest point sampling [not implemented yet]

    For the naive algorithm, max and min values can be specified for value or
    absolute values, if one component exceeds the range for a data point, that
    data point is filtered out.
    """
    import numpy as np
    filterFns = {
        'val_max': lambda data, tol: np.any(data>tol),
        'val_min': lambda data, tol: np.any(data<tol),
        'abs_max': lambda data, tol: np.any(np.abs(data)>tol),
        'abs_min': lambda data, tol: np.any(np.abs(data)<tol)}
    def filter_fn(data, key, tag, tol):
        if tag.startswith('!'):
            return not filterFns[key](data[f"{tag[1:]}_data"], tol)
        else:
            return filterFns[key](data[f"{tag}_data"], tol)
    writer = get_writer(output, format=oformat)
    fnList = [lambda data,k=key,t=tag,tol=tol: filter_fn(data, k, t, float(tol))
              for key, val in kwargs.items() if val
              for tag, tol in map(lambda x: x.split(':'), val.split(','))]
    idx = []
    for fname in filename:
        ds = read(fname, log=log, format=format, emap=emap, units=units)
        for i,data in enumerate(ds):
            if not any([fn(data) for fn in fnList]):
                writer.add(data)
                idx.append(i)
    np.savetxt(f'{output}.idx', idx, fmt='%i')
    writer.finalize()


@click.command(name='merge', context_settings=CONTEXT_SETTINGS, short_help='filter datasets')
@click.argument('filename', nargs=-1)
@click.option('-o', '--output', metavar='', default='dataset')
@click.option('-of', '--oformat', metavar='', default='pinn', help='output format')
def mergeds(filename, output, oformat):
    writer = get_writer(output, format=oformat)
    for fname in filename:
        ds = read(fname, format='pinn')
        for data in ds:
            writer.add(data)
    writer.finalize()

main.add_command(convertds)
main.add_command(splitds)
main.add_command(filterds)
main.add_command(mergeds)
main.add_command(qbc)
main.add_command(version)

if __name__ == '__main__':
    main()
