#!/usr/bin/env python3
# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Wrapper to run source_spec.py from source tree.

:copyright:
    2016-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import sys
import os
import inspect

MIN_PYTHON_VERSION = (3, 6)
MIN_PYTHON_VERSION_STR = '{}.{}'.format(*MIN_PYTHON_VERSION)
PYTHON_VERSION_STR = '{}.{}.{}'.format(*sys.version_info[0:3])
if sys.version_info < MIN_PYTHON_VERSION:
    msg = 'SourceSpec requires Python version >= {}'.format(
        MIN_PYTHON_VERSION_STR)
    msg += ' you are using Python version {}'.format(PYTHON_VERSION_STR)
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
        import obspy #NOQA
        # Try to import tzlocal, which is required at the end of the run
        import tzlocal #NOQA
        from sourcespec.source_spec import main
        main()
    except ImportError as msg:
        mod_name = msg.name
        if mod_name == 'PIL':
            mod_name = 'pillow'
        s = "Error: module '{}' is required by source_spec.".format(mod_name)
        s += " Please install it.\n"
        sys.stderr.write(s)
        sys.exit(1)
