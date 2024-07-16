# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Augment event with velocity info and event name.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2024 Claudio Satriano <satriano@ipgp.fr>,
              Sophie Lambotte <sophie.lambotte@unistra.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
from ..setup import config, ssp_exit
from ..ssp_util import MediumProperties
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _hypo_vel(hypo):
    """
    Compute velocity at hypocenter.

    :param hypo: Hypocenter object
    :type hypo: :class:`sourcespec.ssp_event.Hypocenter`
    """
    medium_properties = MediumProperties(
        hypo.longitude.value_in_deg, hypo.latitude.value_in_deg,
        hypo.depth.value_in_km
    )
    hypo.vp = medium_properties.get(mproperty='vp', where='source')
    hypo.vs = medium_properties.get(mproperty='vs', where='source')
    hypo.rho = medium_properties.get(mproperty='rho', where='source')
    depth_string = medium_properties.to_string(
        'source depth', hypo.depth.value_in_km)
    vp_string = medium_properties.to_string('vp_source', hypo.vp)
    vs_string = medium_properties.to_string('vs_source', hypo.vs)
    rho_string = medium_properties.to_string('rho_source', hypo.rho)
    logger.info(f'{depth_string}, {vp_string}, {vs_string}, {rho_string}')


def augment_event(ssp_event):
    """
    Add velocity info to hypocenter
    and add event name from/to config.options

    The augmented event is stored in config.event

    :param ssp_event: Evento to be augmented
    :type ssp_event: :class:`sourcespec.ssp_event.SSPEvent`
    """
    # add velocity info to hypocenter
    try:
        _hypo_vel(ssp_event.hypocenter)
    except Exception as e:
        logger.error(
            f'Unable to compute velocity at hypocenter: {e}\n')
        ssp_exit(1)
    if config.options.evname is not None:
        # add evname from command line, if any, overriding the one in ssp_event
        ssp_event.name = config.options.evname
    else:
        # add evname from ssp_event, if any, to config file
        config.options.evname = ssp_event.name
    # add event to config file
    config.event = ssp_event
