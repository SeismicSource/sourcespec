# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
An exit function for SourceSpec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import sys
import logging
import signal
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def ssp_exit(retval=0, abort=False):
    """Exit the program."""
    # ssp_exit might have already been called if multiprocessing
    if ssp_exit.SSP_EXIT_CALLED:
        return
    ssp_exit.SSP_EXIT_CALLED = True
    if abort:
        print('\nAborting.')
        if logger is not None:
            logger.debug('source_spec ABORTED')
    elif logger is not None:
        logger.debug('source_spec END')
    logging.shutdown()
    sys.exit(retval)


ssp_exit.SSP_EXIT_CALLED = False


def _sigint_handler(sig, frame):
    """Handle SIGINT signal."""
    # pylint: disable=unused-argument
    ssp_exit(1, abort=True)


signal.signal(signal.SIGINT, _sigint_handler)
