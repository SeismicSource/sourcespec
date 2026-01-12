# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Input functions for SourceSpec.

:copyright:
    2013-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
from .traces import read_traces  # noqa
from .augment_traces import augment_traces  # noqa
from .event_and_picks import read_event_and_picks  # noqa
from .augment_event import augment_event  # noqa
from .station_metadata import read_station_metadata  # noqa
