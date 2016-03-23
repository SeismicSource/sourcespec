# -*- coding: utf-8 -*-
"""setup.py: setuptools control."""
from setuptools import setup

import inspect
import os
import sys

# Import the version string.
path = os.path.join(os.path.abspath(os.path.dirname(inspect.getfile(
    inspect.currentframe()))), 'sourcespec')
sys.path.insert(0, path)
from version import get_git_version


with open('README.md', 'rb') as f:
    long_descr = f.read().decode('utf-8')


setup(
    name='sourcespec',
    packages=['sourcespec'],
    include_package_data=True,
    entry_points={
        'console_scripts': ['source_spec = sourcespec.source_spec:main',
                            'source_model = sourcespec.source_model:main',
                            'source_residuals = '
                            'sourcespec.source_residuals:main']
        },
    version=get_git_version(),
    description='Earthquake source parameters from S-wave '
                'displacement spectra',
    long_description=long_descr,
    author='Claudio Satriano',
    author_email='satriano@gmail.com',
    url='',
    install_requires=['obspy>=1.0.0']
    )
