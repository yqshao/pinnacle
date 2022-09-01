#!/usr/bin/env python3

import os, re
from setuptools import setup, find_packages

with open("tips/__init__.py") as f:
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                      f.read(), re.M)
if match:
    version = match.group(1)
else:
    raise RuntimeError("Unable to find version string.")

setup(
    name="tips",
    version=version,
    description="The Interatomic Potential Suite",
    url="https://github.com/teoroo-CMC/tips",
    author="Yunqi Shao",
    author_email="yunqi_shao@yahoo.com",
    license="BSD",
    packages=find_packages(),
    install_requires=[
        "click>=7.0",
        "numpy>1.3.0",
        "pyyaml>=3.01",
        "ase>=3.22",
        "mock>=4.0"
    ],
    python_requires=">=3.6",
    entry_points={"console_scripts": ["tips=tips.cli.cli:cli"]},
)
