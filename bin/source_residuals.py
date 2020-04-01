#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
Wrapper to run source_residuals.py from source tree.

:copyright:
    2016-2020 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import os
import inspect

# Make sure we use current-dir version over installed one
path = os.path.abspath(os.path.join(os.path.dirname(inspect.getfile(
    inspect.currentframe())), os.pardir))
sys.path.insert(0, path)
from sourcespec.source_residuals import main

if __name__ == '__main__':
    main()
