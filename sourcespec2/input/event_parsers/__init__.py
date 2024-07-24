# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Event and phase picks parsers for SourceSpec.

:copyright:
    2013-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
from .asdf_event import parse_asdf_event_picks  # noqa
from .hypo71 import parse_hypo71_hypocenter, parse_hypo71_picks  # noqa
from .hypo2000 import parse_hypo2000_file  # noqa
from .obspy_catalog import parse_obspy_catalog  # noqa
from .quakeml import parse_qml_event_picks  # noqa
from .source_spec_event import parse_source_spec_event_file  # noqa
from .sac_event import read_event_from_SAC, read_picks_from_SAC  # noqa
