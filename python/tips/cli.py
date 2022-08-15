#!/usr/bin/env python3
import tips
import os
import logging
import click
import click_logging

logger = logging.getLogger('tips')
click_logging.basic_config(logger)

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group()
def main():
    """TIPS CLI - A command line tool for manipulating of AML data"""
    pass


@click.command()
def version():
    click.echo(f"TIPS version: {tips.__version__}")


@click.command(
    name="convert", context_settings=CONTEXT_SETTINGS, short_help="convert datasets"
)
@click_logging.simple_verbosity_option(logger)
@click.argument("filename", nargs=-1)
@click.option("-f", "--fmt", metavar="", default="auto", help="input format")
@click.option("-o", "--output", metavar="", default="dataset", help="output name")
@click.option("-of", "--ofmt", metavar="", default="runner", help="output format")
@click.option("-em", "--emap", metavar="", default=None, help="lammps data for elem")
def convert(filename, fmt, output, ofmt, emap):
    from tips.io import load_ds

    ds = load_ds(filename, fmt=fmt)
    if emap:
        ds.map_elems(emap)
    ds.convert(output, fmt=ofmt)


main.add_command(convert)
main.add_command(version)

if __name__ == "__main__":
    main()
