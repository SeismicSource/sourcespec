# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Setup functions for sourcespec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import sys
import os
import platform
import shutil
import logging
import signal
import contextlib
from datetime import datetime
from . import __version__, __banner__
from .config import config, library_versions

# define ipshell(), if possible
# note: ANSI colors do not work on Windows standard terminal
if sys.stdout.isatty() and sys.platform != 'win32':
    try:
        from IPython.terminal.embed import InteractiveShellEmbed
        IPSHELL = InteractiveShellEmbed.instance()
    except ImportError:
        IPSHELL = None
else:
    IPSHELL = None

# global variables
OLDLOGFILE = None
LOGGER = None
SSP_EXIT_CALLED = False


def save_config():
    """Save config file to output dir."""
    # Actually, it renames the file already existing.
    src = os.path.join(config.options.outdir, 'source_spec.conf')
    evid = config.event.event_id
    dst = os.path.join(config.options.outdir, f'{evid}.ssp.conf')
    # On Windows, dst file must not exist
    with contextlib.suppress(Exception):
        os.remove(dst)
    os.rename(src, dst)


def move_outdir():
    """Move outdir to a new dir named from evid (and optional run_id)."""
    try:
        evid = config.event.event_id
    except Exception:
        return
    src = config.options.outdir
    run_id = config.options.run_id
    run_id_subdir = config.options.run_id_subdir
    dst = os.path.split(src)[0]
    dst = os.path.join(dst, str(evid))
    if run_id and run_id_subdir:
        dst = os.path.join(dst, str(run_id))
    elif run_id:
        dst += f'_{run_id}'
    # Create destination
    if not os.path.exists(dst):
        os.makedirs(dst)
    # Copy all files into destination
    file_names = os.listdir(src)
    for file_name in file_names:
        shutil.copyfile(
            os.path.join(src, file_name),
            os.path.join(dst, file_name)
        )
    # Old outdir cannot be removed yet, because the log file is still opened
    config.options.oldoutdir = src
    config.options.outdir = dst


def remove_old_outdir():
    """Try to remove the old outdir."""
    try:
        oldoutdir = config.options.oldoutdir
        shutil.rmtree(oldoutdir)
    except Exception:
        return


def _color_handler_emit(fn):
    """
    Add color-coding to the logging handler emitter.

    Source: https://stackoverflow.com/a/20707569/2021880
    """
    def new(*args):
        levelno = args[0].levelno
        if levelno >= logging.CRITICAL:
            color = '\x1b[31;1m'  # red
        elif levelno >= logging.ERROR:
            color = '\x1b[31;1m'  # red
        elif levelno >= logging.WARNING:
            color = '\x1b[33;1m'  # yellow
        elif levelno >= logging.INFO:
            # color = '\x1b[32;1m'  # green
            color = '\x1b[0m'  # no color
        elif levelno >= logging.DEBUG:
            color = '\x1b[35;1m'  # purple
        else:
            color = '\x1b[0m'  # no color
        # Color-code the message
        args[0].msg = f'{color}{args[0].msg}\x1b[0m'
        return fn(*args)
    return new


def setup_logging(basename=None, progname='source_spec'):
    """
    Set up the logging infrastructure.

    This function is typically called twice: the first time without basename
    and a second time with a basename (typically the eventid).
    When called the second time, the previous logfile is renamed using the
    given basename.

    :param basename: The basename for the logfile (default: None).
    :type basename: str
    :param progname: The program name to be used in the log file (default:
        'source_spec')
    :type progname: str
    """
    # pylint: disable=global-statement
    global OLDLOGFILE
    global LOGGER
    # Create outdir
    if not os.path.exists(config.options.outdir):
        os.makedirs(config.options.outdir)

    if basename:
        logfile = os.path.join(config.options.outdir, f'{basename}.ssp.log')
    else:
        datestring = datetime.now().strftime('%Y%m%d_%H%M%S')
        logfile = os.path.join(config.options.outdir, f'{datestring}.ssp.log')

    logger_root = logging.getLogger()
    if OLDLOGFILE:
        hdlrs = logger_root.handlers[:]
        for hdlr in hdlrs:
            hdlr.flush()
            hdlr.close()
            logger_root.removeHandler(hdlr)
        # Copy old logfile to new
        shutil.copyfile(OLDLOGFILE, logfile)
        # Remove old logfile from old and new dir
        os.remove(OLDLOGFILE)
        OLDLOGFILE = os.path.join(
            config.options.outdir, os.path.basename(OLDLOGFILE))
        os.remove(OLDLOGFILE)
        filemode = 'a'
    else:
        filemode = 'w'

    # captureWarnings is not supported in old versions of python
    with contextlib.suppress(Exception):
        logging.captureWarnings(True)
    logger_root.setLevel(logging.DEBUG)
    filehand = logging.FileHandler(filename=logfile, mode=filemode)
    filehand.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s %(name)-20s '
                                  '%(levelname)-8s %(message)s')
    filehand.setFormatter(formatter)
    logger_root.addHandler(filehand)

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)

    # Add logger color coding on all platforms but win32
    if sys.platform != 'win32' and sys.stdout.isatty():
        console.emit = _color_handler_emit(console.emit)
    logger_root.addHandler(console)

    LOGGER = logging.getLogger(progname)

    # Only write these debug infos for a new logfile
    if not OLDLOGFILE:
        _log_debug_information()
    OLDLOGFILE = logfile

    # See if there are warnings to deliver
    for _ in range(len(config.warnings)):
        msg = config.warnings.pop(0)
        LOGGER.warning(msg)


def _log_debug_information():
    banner = f'\n{__banner__}\nThis is SourceSpec v{__version__}.\n'
    LOGGER.info(banner)
    LOGGER.debug('source_spec START')
    LOGGER.debug(f'SourceSpec version: {__version__}')
    uname = platform.uname()
    uname_str = f'{uname[0]} {uname[2]} {uname[4]}'
    LOGGER.debug(f'Platform: {uname_str}')
    lv = library_versions
    LOGGER.debug(f'Python version: {lv.PYTHON_VERSION_STR}')
    LOGGER.debug(f'ObsPy version: {lv.OBSPY_VERSION_STR}')
    LOGGER.debug(f'NumPy version: {lv.NUMPY_VERSION_STR}')
    LOGGER.debug(f'SciPy version: {lv.SCIPY_VERSION_STR}')
    LOGGER.debug(f'Matplotlib version: {lv.MATPLOTLIB_VERSION_STR}')
    LOGGER.debug(f'PyProj version: {lv.PYPROJ_VERSION_STR}')
    if lv.CARTOPY_VERSION_STR is not None:
        LOGGER.debug(f'Cartopy version: {lv.CARTOPY_VERSION_STR}')
    LOGGER.debug('Running arguments:')
    LOGGER.debug(' '.join(sys.argv))


def ssp_exit(retval=0, abort=False):
    """Exit the program."""
    # ssp_exit might have already been called if multiprocessing
    global SSP_EXIT_CALLED  # pylint: disable=global-statement
    if SSP_EXIT_CALLED:
        return
    SSP_EXIT_CALLED = True
    if abort:
        print('\nAborting.')
        if LOGGER is not None:
            LOGGER.debug('source_spec ABORTED')
    elif LOGGER is not None:
        LOGGER.debug('source_spec END')
    logging.shutdown()
    sys.exit(retval)


def sigint_handler(sig, frame):
    """Handle SIGINT signal."""
    # pylint: disable=unused-argument
    ssp_exit(1, abort=True)


signal.signal(signal.SIGINT, sigint_handler)
