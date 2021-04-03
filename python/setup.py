#!/usr/bin/env python3

import os, re
from setuptools import setup, find_packages

with open('tips/__init__.py') as f:
    version = re.search("__version__ = '(.*)'", f.read()).group(1)

setup(name='tips',
      version=version,
      description='Pair interaction neural network',
      url='https://github.com/teoroo-CMC/pinn',
      author='Yunqi Shao',
      author_email='yunqi_shao@yahoo.com',
      license='BSD',
      packages=find_packages(),
      install_requires=['click>=7.0',
                        'numpy>1.3.0',
                        'ase>=3.19.0',
                        'pyyaml>=3.01'],
      python_requires='>=3.6',
      entry_points={'console_scripts':
                    ['tips=tips.cli:main']}
)
