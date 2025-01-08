# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
An exit function for SourceSpec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import sys
import logging
import signal
from . import config
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def ssp_exit(status=0, abort=False, extra_message=None):
    """
    Exit the program.

    If the status is an integer, it is used as the exit status.
    If the status is a string, it is used as an error message and the exit
    status is set to 1 (failure).

    :param status: Exit status
    :type status: int or str
    :param abort: If True, print 'Aborting.' before exiting
    :type abort: bool
    :param extra_message: Extra message to print before exiting.
        It will be printed to stderr.
    :type extra_message: str
    """
    # ssp_exit might have already been called if multiprocessing
    if ssp_exit.SSP_EXIT_CALLED:
        return
    ssp_exit.SSP_EXIT_CALLED = True
    if isinstance(status, str):
        err_message = status
        logger.error(err_message)
        retval = 1
    else:
        err_message = None
        retval = status
    if abort:
        print('\nAborting.')
        if logger is not None:
            logger.debug('source_spec ABORTED')
    elif logger is not None:
        logger.debug('source_spec END')
    logging.shutdown()
    if not config.running_from_command_line and retval:
        raise RuntimeError(err_message)
    if extra_message:
        sys.stderr.write(extra_message)
    sys.exit(retval)


ssp_exit.SSP_EXIT_CALLED = False


def _sigint_handler(sig, frame):
    """Handle SIGINT signal."""
    # pylint: disable=unused-argument
    ssp_exit(1, abort=True)


signal.signal(signal.SIGINT, _sigint_handler)
