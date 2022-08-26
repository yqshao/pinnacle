#!/usr/bin/env python
import click
from .main import convert, subsample, version
from .utils import mkcp2kinp

@click.group()
def cli():
    """TIPS CLI - A command line tool for manipulating of AML data"""
    pass
cli.add_command(convert)
cli.add_command(subsample)
cli.add_command(version)

@cli.group()
def utils():
    pass
utils.add_command(mkcp2kinp)


if __name__ == "__main__":
    cli()
