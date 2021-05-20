# -*- coding: utf-8 -*-
"""setup.py: setuptools control."""
from setuptools import setup
import versioneer

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
    long_description='Earthquake source parameters from S-wave '
                'displacement spectra',
    long_description_content_type='text/markdown',
    author='Claudio Satriano',
    author_email='satriano@ipgp.fr',
    url='http://www.ipgp.fr/~satriano',
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
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Physics'],
    install_requires=['obspy>=1.1.0', 'scipy>=0.17', 'matplotlib>=2.2', 'pytz']
    )
