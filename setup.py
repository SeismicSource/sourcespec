# -*- coding: utf-8 -*-
"""setup.py: setuptools control."""

from setuptools import setup

version = 0.9

with open('README.md', 'rb') as f:
    long_descr = f.read().decode('utf-8')


setup(
    name='source_spec',
    packages=['source_spec'],
    entry_points={
        'console_scripts': ['source_spec = source_spec.source_spec:main',
                            'source_model = source_spec.source_model:main',
                            'source_residuals = '
                            'source_spec.source_residuals:main']
        },
    version=version,
    description='Modelling S-wave displacement spectra and '
                'inverting source parameters.',
    long_description=long_descr,
    author='Claudio Satriano',
    author_email='satriano@gmail.com',
    url='',
    )
