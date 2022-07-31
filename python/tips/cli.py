#!/usr/bin/env python3
import tips
import click

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
@click.argument("filename", nargs=-1)
@click.option("-f", "--fmt", metavar="", default="auto", help="input format")
@click.option("-o", "--output", metavar="", default="dataset")
@click.option("-of", "--ofmt", metavar="", default="runner", help="output format")
def convert(filename, fmt, output, ofmt):
    from tips.io import load_ds
    ds = load_ds(filename, fmt=fmt)
    ds.convert(output, fmt=ofmt)


main.add_command(convert)
main.add_command(version)

if __name__ == "__main__":
    main()
