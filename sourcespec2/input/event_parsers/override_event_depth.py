# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Override event depth from command line argument.

:copyright:
    2012-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def override_event_depth(ssp_event, depth_in_km):
    """
    Override event depth if specified in command line arguments.

    :param ssp_event: SSPEvent object
    :param depth_in_km: depth value in kilometers
    """
    if depth_in_km is not None:
        ssp_event.hypocenter.depth.value = depth_in_km
        ssp_event.hypocenter.depth.units = 'km'
        logger.info(f'Overriding event depth to {depth_in_km} km')
