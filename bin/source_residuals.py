#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Wrapper to run source_residuals.py from source tree."""
import sys
import os
import inspect

path = os.path.abspath(os.path.join(os.path.dirname(inspect.getfile(
    inspect.currentframe())), os.pardir, 'sourcespec'))
sys.path.insert(0, path)
from source_residuals import main

if __name__ == '__main__':
    main()
