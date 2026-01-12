# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Configuration classes and functions for SourceSpec

:copyright:
    2013-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
from .config import config  # noqa
from .configure_cli import configure_cli  # noqa
from .logging import setup_logging  # noqa
from .exit import ssp_exit  # noqa
from .outdir import (  # noqa
    get_outdir_path, save_config, move_outdir, remove_old_outdir
)
