#!/usr/bin/env python3
# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Wrapper to run source_spec.py from source tree.

:copyright:
    2016-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import sys
import os
import inspect

MIN_PYTHON_VERSION = (3, 6)
# pylint: disable=consider-using-f-string
MIN_PYTHON_VERSION_STR = '{}.{}'.format(*MIN_PYTHON_VERSION)
PYTHON_VERSION_STR = '{}.{}.{}'.format(*sys.version_info[:3])
if sys.version_info < MIN_PYTHON_VERSION:
    msg = f'SourceSpec requires Python version >= {MIN_PYTHON_VERSION_STR}'
    msg += f' you are using Python version {PYTHON_VERSION_STR}'
    print(msg, file=sys.stderr)
    sys.exit(1)

if __name__ == '__main__':
    try:
        # Make sure we use current-dir version over installed one
        path = os.path.abspath(os.path.join(os.path.dirname(inspect.getfile(
            inspect.currentframe())), os.pardir))
        sys.path.insert(0, path)
        # Try to import obspy, which requires most of the
        # source_spec dependencies
        import obspy  # NOQA  pylint: disable=unused-import
        # Try to import tzlocal, which is required at the end of the run
        import tzlocal  # NOQA  pylint: disable=unused-import
        from sourcespec.source_spec import main
        main()
    except ImportError as msg:
        MOD_NAME = msg.name
        if MOD_NAME == 'PIL':
            MOD_NAME = 'pillow'
        sys.stderr.write(
            f"Error: module '{MOD_NAME}' is required by source_spec. "
            "Please install it.\n"
        )
        sys.exit(1)
