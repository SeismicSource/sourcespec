# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Retrieve and check library versions for SourceSpec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import sys
import os
import contextlib
import warnings


# ---- Helper functions for cartopy feature download ----
def _cartopy_download_gshhs():
    """
    Download GSHHS data for cartopy.
    """
    # pylint: disable=import-outside-toplevel
    from cartopy.io.shapereader import GSHHSShpDownloader as Downloader
    from cartopy import config as cartopy_config
    from pathlib import Path
    gshhs_downloader = Downloader.from_config(('shapefiles', 'gshhs'))
    format_dict = {'config': cartopy_config, 'scale': 'f', 'level': 1}
    target_path = gshhs_downloader.target_path(format_dict)
    if not os.path.exists(target_path):
        sys.stdout.write(
            'Downloading GSHHS data for cartopy.\n'
            'This is needed only the first time you use cartopy and may take '
            'a while...\n')
        sys.stdout.flush()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            try:
                path = gshhs_downloader.path(format_dict)
            except Exception:
                sys.stderr.write(
                    '\nUnable to download data. '
                    'Check your internet connection.\n')
                sys.exit(1)
        sys.stdout.write(f'Done! Data cached to {Path(path).parents[1]}\n\n')


def _cartopy_download_borders():
    """
    Download borders data for cartopy.

    Inspired from
    https://github.com/SciTools/cartopy/blob/main/tools/cartopy_feature_download.py
    """
    # pylint: disable=import-outside-toplevel
    from cartopy.io import Downloader
    from cartopy import config as cartopy_config
    from pathlib import Path
    category = 'cultural'
    name = 'admin_0_boundary_lines_land'
    scales = ('10m', '50m')
    for scale in scales:
        downloader = Downloader.from_config((
            'shapefiles', 'natural_earth', category, name, scale))
        format_dict = {
            'config': cartopy_config, 'category': category, 'name': name,
            'resolution': scale}
        target_path = downloader.target_path(format_dict)
        if not os.path.exists(target_path):
            sys.stdout.write(
                'Downloading border data for cartopy...\n')
            sys.stdout.flush()
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                try:
                    path = downloader.path(format_dict)
                except Exception:
                    sys.stderr.write(
                        '\nUnable to download data. '
                        'Check your internet connection.\n')
                    sys.exit(1)
            sys.stdout.write(
                f'Done! Data cached to {Path(path).parents[0]}\n\n')


class _LibraryVersions():
    """
    Check library versions and download data for cartopy.

    This is a private class, only the instance `library_versions` is
    accessible from the outside.
    """
    def __init__(self):
        self.MAX_NUMPY_VERSION = (2, 0, 0)
        self.MAX_MATPLOTLIB_VERSION = (3, 10, 0)
        self.MIN_OBSPY_VERSION = (1, 2, 0)
        self.MIN_CARTOPY_VERSION = (0, 21, 0)
        self.MIN_NLLGRID_VERSION = (1, 4, 2)
        self.OBSPY_VERSION = None
        self.OBSPY_VERSION_STR = None
        self.NUMPY_VERSION_STR = None
        self.SCIPY_VERSION_STR = None
        self.MATPLOTLIB_VERSION_STR = None
        self.CARTOPY_VERSION_STR = None
        self.PYTHON_VERSION_STR = None
        # check base library versions and obspy version
        # the other checks are called on demand
        self.check_base_library_versions()
        self.check_obspy_version()

    def check_base_library_versions(self):
        """
        Check base library versions.
        """
        # pylint: disable=import-outside-toplevel
        self.PYTHON_VERSION_STR = '.'.join(map(str, sys.version_info[:3]))
        import numpy
        self.NUMPY_VERSION_STR = numpy.__version__
        numpy_version = tuple(map(int, self.NUMPY_VERSION_STR.split('.')[:3]))
        if numpy_version >= self.MAX_NUMPY_VERSION:
            max_numpy_version_str = '.'.join(map(str, self.MAX_NUMPY_VERSION))
            raise ImportError(
                f'ERROR: Numpy >= {max_numpy_version_str} '
                'is not yet supported. Please use a less recent version. '
                f'You have version: {self.NUMPY_VERSION_STR}'
            )
        import scipy
        self.SCIPY_VERSION_STR = scipy.__version__
        import matplotlib
        self.MATPLOTLIB_VERSION_STR = matplotlib.__version__
        matplotlib_version = self.MATPLOTLIB_VERSION_STR.split('.')[:3]
        matplotlib_version = tuple(map(int, matplotlib_version))
        if matplotlib_version >= self.MAX_MATPLOTLIB_VERSION:
            max_matplotlib_version_str =\
                '.'.join(map(str, self.MAX_MATPLOTLIB_VERSION))
            raise ImportError(
                f'ERROR: Matplotlib >= {max_matplotlib_version_str} '
                'is not yet supported. Please use a less recent version'
                f' You have version: {self.MATPLOTLIB_VERSION_STR}'
            )

    def check_obspy_version(self):
        """
        Check ObsPy version.
        """
        # pylint: disable=import-outside-toplevel
        import obspy
        self.OBSPY_VERSION_STR = obspy.__version__
        obspy_version = self.OBSPY_VERSION_STR.split('.')[:3]
        # special case for "rc" versions:
        obspy_version[2] = obspy_version[2].split('rc')[0]
        obspy_version = tuple(map(int, obspy_version))
        with contextlib.suppress(IndexError):
            # add half version number for development versions
            # check if there is a fourth field in version string:
            obspy_version = obspy_version[:2] + (obspy_version[2] + 0.5,)
        if obspy_version < self.MIN_OBSPY_VERSION:
            min_obspy_version_str = '.'.join(map(str, self.MIN_OBSPY_VERSION))
            raise ImportError(
                f'ERROR: ObsPy >= {min_obspy_version_str} is required. '
                f'You have version: {self.OBSPY_VERSION_STR}'
            )

    def check_cartopy_version(self):
        """
        Check cartopy version and download data if needed.
        """
        try:
            cartopy_ver = None
            import cartopy  # NOQA pylint: disable=import-outside-toplevel
            self.CARTOPY_VERSION_STR = cartopy.__version__
            cartopy_ver = tuple(map(int, cartopy.__version__.split('.')[:3]))
            if cartopy_ver < self.MIN_CARTOPY_VERSION:
                raise ImportError
            _cartopy_download_gshhs()
            _cartopy_download_borders()
        except ImportError as e:
            cartopy_min_ver_str = '.'.join(map(str, self.MIN_CARTOPY_VERSION))
            msg = (
                f'\nPlease install cartopy >= {cartopy_min_ver_str} '
                'to plot maps.\nHow to install: '
                'https://scitools.org.uk/cartopy/docs/latest/installing.html\n'
                '\nAlternatively, set "plot_station_map" to "False" '
                'in config file.\n'
            )
            if cartopy_ver is not None:
                msg += (
                    f'Installed cartopy version: {self.CARTOPY_VERSION_STR}.\n'
                )
            raise ImportError(msg) from e

    def check_rasterio_version(self):
        """
        Check rasterio version.
        """
        # pylint: disable=import-outside-toplevel
        try:
            import rasterio  # noqa pylint: disable=unused-import
        except ImportError as e:
            msg = (
                '\nPlease install rasterio to plot GeoTIFF files.\n'
                'How to install: https://rasterio.readthedocs.io/en/stable/\n'
            )
            raise ImportError(msg) from e

    def check_pyproj_version(self):
        """
        Check pyproj version.
        """
        # pylint: disable=import-outside-toplevel
        try:
            import pyproj # noqa pylint: disable=unused-import
        except ImportError as e:
            msg = '\nPlease install pyproj to plot maps.\n'
            raise ImportError(msg) from e

    def check_nllgrid_version(self):
        """
        Check nllgrid version.
        """
        # pylint: disable=import-outside-toplevel
        try:
            nllgrid_ver = None
            import nllgrid  # NOQA
            nllgrid_ver_str = nllgrid.__version__.split('+')[0]
            nllgrid_ver = tuple(map(int, nllgrid_ver_str.split('.')))
            # nllgrid versions are sometimes X.Y, other times X.Y.Z
            while len(nllgrid_ver) < 3:
                nllgrid_ver = (*nllgrid_ver, 0)
            if nllgrid_ver < self.MIN_NLLGRID_VERSION:
                raise ImportError
        except ImportError as e:
            nllgrid_min_ver_str = '.'.join(map(str, self.MIN_NLLGRID_VERSION))
            msg = (
                f'\nPlease install nllgrid >= {nllgrid_min_ver_str} to use '
                'NonLinLoc grids.\n'
                'How to install: https://github.com/claudiodsf/nllgrid\n'
            )
            if nllgrid_ver is not None:
                msg += f'Installed nllgrid version: {nllgrid.__version__}\n'
            raise ImportError(msg) from e


# class instance exposed to the outside
try:
    library_versions = _LibraryVersions()
except ImportError as _error:
    sys.exit(_error)
