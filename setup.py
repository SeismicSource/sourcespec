# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
A minimal setup script for SourceSpec.

This script ensures compatibility with versioneer
and dynamically generates a README for correctly displaying images on PyPI.

All the remaining configuration is in pyproject.toml.
"""
from setuptools import setup
import versioneer

# Dynamically generate the README text for PyPI, replacing local image paths
# with the corresponding URLs on the CDN.
revision = versioneer.get_versions()['full-revisionid']
cdn_baseurl =\
    f'https://cdn.jsdelivr.net/gh/SeismicSource/sourcespec@{revision}'
with open('README.md', 'rb') as f:
    long_description = f.read().decode('utf-8').replace(
        'imgs/SourceSpec_logo.svg',
        f'{cdn_baseurl}/imgs/SourceSpec_logo.svg'
    ).replace(
        'imgs/example_trace.svg',
        f'{cdn_baseurl}/imgs/example_trace.svg'
    ).replace(
        'imgs/example_spectrum.svg',
        f'{cdn_baseurl}/imgs/example_spectrum.svg'
    ).replace(
        '(CHANGELOG.md)',
        f'({cdn_baseurl}/CHANGELOG.md)'
    )

setup(
    long_description=long_description,
    long_description_content_type='text/markdown',
    cmdclass=versioneer.get_cmdclass()
)
