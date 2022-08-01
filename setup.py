# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""setup.py: setuptools control."""
from setuptools import setup
import versioneer
import sys

MIN_PYTHON_VERSION = (3, 6)
MIN_PYTHON_VERSION_STR = '{}.{}'.format(*MIN_PYTHON_VERSION)
PYTHON_VERSION_STR = '{}.{}.{}'.format(*sys.version_info[0:3])
if sys.version_info < MIN_PYTHON_VERSION:
    msg = 'SourceSpec requires Python version >= {}'.format(
        MIN_PYTHON_VERSION_STR)
    msg += ' you are using Python version {}'.format(PYTHON_VERSION_STR)
    print(msg, file=sys.stderr)
    sys.exit(1)

revision = versioneer.get_versions()['full-revisionid']
cdn_baseurl = 'https://cdn.jsdelivr.net/gh/SeismicSource/sourcespec@{}'\
    .format(revision)
with open('README.md', 'rb') as f:
    long_descr = f.read().decode('utf-8').replace(
        'imgs/SourceSpec_logo.svg',
        '{}/imgs/SourceSpec_logo.svg'.format(cdn_baseurl)
    ).replace(
        'imgs/example_trace.svg',
        '{}/imgs/example_trace.svg'.format(cdn_baseurl)
    ).replace(
        'imgs/example_spectrum.svg',
        '{}/imgs/example_spectrum.svg'.format(cdn_baseurl)
    ).replace(
        '(CHANGELOG.md)',
        '({}/CHANGELOG.md)'.format(cdn_baseurl)
    )

project_urls = {
    'Homepage': 'https://sourcespec.seismicsource.org',
    'Source': 'https://github.com/SeismicSource/sourcespec',
    'Documentation': 'https://sourcespec.readthedocs.io'
}

setup(
    name='sourcespec',
    packages=['sourcespec', 'sourcespec.configobj', 'sourcespec.adjustText'],
    include_package_data=True,
    entry_points={
        'console_scripts': ['source_spec = sourcespec.source_spec:main',
                            'source_model = sourcespec.source_model:main',
                            'source_residuals = '
                            'sourcespec.source_residuals:main']
        },
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Earthquake source parameters from S-wave '
                'displacement spectra',
    long_description=long_descr,
    long_description_content_type='text/markdown',
    author='Claudio Satriano',
    author_email='satriano@ipgp.fr',
    project_urls=project_urls,
    license='CeCILL Free Software License Agreement, Version 2.1',
    platforms='OS Independent',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: CEA CNRS Inria Logiciel Libre '
            'License, version 2.1 (CeCILL-2.1)',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics'],
    python_requires='>={}'.format(MIN_PYTHON_VERSION_STR),
    install_requires=[
        'numpy>=1.10',
        'scipy>=0.17',
        'matplotlib>=2.2',
        'pillow>=4.0.0',
        'obspy>=1.2.0',
        'pyproj',
        'tzlocal']
    )
