#!/usr/bin/env python3
# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Wrapper to run source_model.py from source tree.

:copyright:
    2016-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import sys
import os
import inspect

if __name__ == '__main__':
    try:
        # Make sure we use current-dir version over installed one
        path = os.path.abspath(os.path.join(os.path.dirname(inspect.getfile(
            inspect.currentframe())), os.pardir))
        sys.path.insert(0, path)
        from sourcespec.source_model import main
        main()
    except ImportError as msg:
        mod_name = msg.name
        s = "Error: module '{}' is required by source_model.".format(mod_name)
        s += " Please install it.\n"
        sys.stderr.write(s)
        sys.exit(1)
