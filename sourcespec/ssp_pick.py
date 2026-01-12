# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
SourceSpec pick class.

:copyright:
    2023-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""


class SSPPick():
    """A pick object."""

    def __init__(self):
        self.station = None
        self.flag = None
        self.phase = None
        self.polarity = None
        self.quality = None
        self.time = None

    def __str__(self):
        return (
            f'station: {self.station}, flag: {self.flag}, '
            f'phase: {self.phase}, polarity: {self.polarity}, '
            f'quality: {self.quality}, time: {self.time}'
        )
