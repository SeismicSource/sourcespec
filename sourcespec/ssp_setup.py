# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Setup functions for sourcespec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2022 Claudio Satriano <satriano@ipgp.fr>
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
import uuid
from datetime import datetime
from collections import defaultdict
from sourcespec.configobj import ConfigObj
from sourcespec.configobj.validate import Validator
from sourcespec.config import Config
from sourcespec._version import get_versions

# define ipshell(), if possible
# note: ANSI colors do not work on Windows standard terminal
if sys.stdout.isatty() and sys.platform != 'win32':
    try:
        from IPython.terminal.embed import InteractiveShellEmbed
        ipshell = InteractiveShellEmbed()
    except ImportError:
        ipshell = None
else:
    ipshell = None

# global variables
OS = os.name
oldlogfile = None
logger = None
ssp_exit_called = False
OBSPY_VERSION = None
OBSPY_VERSION_STR = None
NUMPY_VERSION_STR = None
SCIPY_VERSION_STR = None
MATPLOTLIB_VERSION_STR = None
CARTOPY_VERSION_STR = None
PYTHON_VERSION_STR = None


def _check_obspy_version():
    global OBSPY_VERSION, OBSPY_VERSION_STR
    # check ObsPy version
    import obspy
    MIN_OBSPY_VERSION = (1, 2, 0)
    OBSPY_VERSION_STR = obspy.__version__
    OBSPY_VERSION = OBSPY_VERSION_STR.split('.')[:3]
    # special case for "rc" versions:
    OBSPY_VERSION[2] = OBSPY_VERSION[2].split('rc')[0]
    OBSPY_VERSION = tuple(map(int, OBSPY_VERSION))
    try:
        # add half version number for development versions
        # check if there is a fourth field in version string:
        OBSPY_VERSION_STR.split('.')[3]
        OBSPY_VERSION = OBSPY_VERSION[:2] + (OBSPY_VERSION[2] + 0.5,)
    except IndexError:
        pass
    if OBSPY_VERSION < MIN_OBSPY_VERSION:
        msg = 'ERROR: ObsPy >= %s.%s.%s is required.' % MIN_OBSPY_VERSION
        msg += ' You have version: %s\n' % OBSPY_VERSION_STR
        sys.stderr.write(msg)
        sys.exit(1)


def _check_cartopy_version():
    cartopy_min_ver = (0, 18, 0)
    try:
        cartopy_ver = None
        import cartopy  #NOQA
        global CARTOPY_VERSION_STR
        CARTOPY_VERSION_STR = cartopy.__version__
        cartopy_ver = tuple(map(int, cartopy.__version__.split('.')[:3]))
        if cartopy_ver < cartopy_min_ver:
            raise ImportError
    except ImportError:
        msg = (
            '\nPlease install cartopy >= {}.{}.{} to plot maps.\n'
            'How to install: '
            'https://scitools.org.uk/cartopy/docs/latest/installing.html\n\n'
            'Alternatively, set "plot_station_map" and "html_report" to '
            '"False" in config file.\n'
            .format(*cartopy_min_ver)
        )
        if cartopy_ver is not None:
            msg += (
                'Installed cartopy version: {}.\n'.format(CARTOPY_VERSION_STR)
            )
        raise ImportError(msg)


def _check_pyproj_version():
    try:
        import pyproj  #NOQA
    except ImportError:
        msg = '\nPlease install pyproj to plot maps.\n'
        raise ImportError(msg)


def _check_nllgrid_version():
    nllgrid_min_ver = (1, 4, 2)
    try:
        nllgrid_ver = None
        import nllgrid  #NOQA
        nllgrid_ver_str = nllgrid.__version__.split('+')[0]
        nllgrid_ver = tuple(map(int, nllgrid_ver_str.split('.')))
        # nllgrid versions are sometimes X.Y, other times X.Y.Z
        while len(nllgrid_ver) < 3:
            nllgrid_ver = (*nllgrid_ver, 0)
        if nllgrid_ver < nllgrid_min_ver:
            raise ImportError
    except ImportError:
        msg = (
            '\nPlease install nllgrid >= {}.{}.{} to use NonLinLoc grids.\n'
            'How to install: https://github.com/claudiodsf/nllgrid\n'
            .format(*nllgrid_min_ver)
        )
        if nllgrid_ver is not None:
            msg += (
                'Installed nllgrid version: {}\n'.format(nllgrid.__version__)
            )
        raise ImportError(msg)


def _check_library_versions():
    global PYTHON_VERSION_STR
    global NUMPY_VERSION_STR
    global SCIPY_VERSION_STR
    global MATPLOTLIB_VERSION_STR
    PYTHON_VERSION_STR = '{}.{}.{}'.format(*sys.version_info[0:3])
    import numpy
    NUMPY_VERSION_STR = numpy.__version__
    import scipy
    SCIPY_VERSION_STR = scipy.__version__
    import matplotlib
    MATPLOTLIB_VERSION_STR = matplotlib.__version__
    MATPLOTLIB_VERSION = MATPLOTLIB_VERSION_STR.split('.')[:3]
    MATPLOTLIB_VERSION = tuple(map(int, MATPLOTLIB_VERSION))
    MAX_MATPLOTLIB_VERSION = (3, 9, 0)
    if MATPLOTLIB_VERSION >= MAX_MATPLOTLIB_VERSION:
        msg = 'ERROR: Matplotlib >= %s.%s.%s ' % MAX_MATPLOTLIB_VERSION
        msg += 'is not yet supported. Please use a less recent version'
        msg += ' You have version: %s\n' % MATPLOTLIB_VERSION_STR
        sys.stderr.write(msg)
        sys.exit(1)


def _read_config(config_file, configspec=None):
    kwargs = dict(
        configspec=configspec, file_error=True, default_encoding='utf8')
    if configspec is None:
        kwargs.update(
            dict(interpolation=False, list_values=False, _inspec=True))
    try:
        config_obj = ConfigObj(config_file, **kwargs)
    except IOError as err:
        sys.stderr.write('{}\n'.format(err))
        sys.exit(1)
    except Exception as err:
        sys.stderr.write('Unable to read "{}": {}\n'.format(config_file, err))
        sys.exit(1)
    return config_obj


def _parse_configspec():
    configspec_file = os.path.join(
        os.path.dirname(__file__), 'configspec.conf')
    configspec = _read_config(configspec_file)
    return configspec


def _write_sample_config(configspec, progname):
    c = ConfigObj(configspec=configspec, default_encoding='utf8')
    val = Validator()
    c.validate(val)
    c.defaults = []
    c.initial_comment = configspec.initial_comment
    c.comments = configspec.comments
    c.final_comment = configspec.final_comment
    configfile = progname + '.conf'
    write_file = True
    if os.path.exists(configfile):
        ans = input(
            '{} already exists. Do you want to overwrite it? [y/N] '.format(
                configfile))
        if ans in ['y', 'Y']:
            write_file = True
        else:
            write_file = False
    if write_file:
        with open(configfile, 'wb') as fp:
            c.write(fp)
        print('Sample config file written to: ' + configfile)


def _update_config_file(config_file, configspec):
    config_obj = _read_config(config_file, configspec)
    val = Validator()
    config_obj.validate(val)
    mod_time = datetime.fromtimestamp(os.path.getmtime(config_file))
    mod_time_str = mod_time.strftime('%Y%m%d_%H%M%S')
    config_file_old = '{}.{}'.format(config_file, mod_time_str)
    ans = input(
        'Ok to update {}? [y/N]\n(Old file will be saved as {}) '.format(
            config_file, config_file_old))
    if ans not in ['y', 'Y']:
        sys.exit(0)
    config_new = ConfigObj(configspec=configspec, default_encoding='utf8')
    config_new = _read_config(None, configspec)
    config_new.validate(val)
    config_new.defaults = []
    config_new.comments = configspec.comments
    config_new.initial_comment = config_obj.initial_comment
    config_new.final_comment = configspec.final_comment
    for k, v in config_obj.items():
        if k not in config_new:
            continue
        # Fix for force_list(default=None)
        if v == ['None', ]:
            v = None
        config_new[k] = v
    migrate_options = {
        's_win_length': 'win_length',
        'traceids': 'traceid_mapping_file',
        'ignore_stations': 'ignore_traceids',
        'use_stations': 'use_traceids',
        'dataless': 'station_metadata',
        'clip_nmax': 'clip_max_percent',
        'PLOT_SHOW': 'plot_show',
        'PLOT_SAVE': 'plot_save',
        'PLOT_SAVE_FORMAT': 'plot_save_format',
        'vp': 'vp_source',
        'vs': 'vs_source',
        'pre_p_time': 'noise_pre_time',
        'pre_s_time': 'signal_pre_time',
        'rps_from_focal_mechanism': 'rp_from_focal_mechanism',
        'paz': 'station_metadata',
    }
    for old_opt, new_opt in migrate_options.items():
        if old_opt in config_obj and config_obj[old_opt] != 'None':
            config_new[new_opt] = config_obj[old_opt]
    shutil.copyfile(config_file, config_file_old)
    with open(config_file, 'wb') as fp:
        config_new.write(fp)
        print('{}: updated'.format(config_file))


def _write_config(config_obj, progname, outdir):
    if progname != 'source_spec':
        return
    configfile = progname + '.conf'
    configfile = os.path.join(outdir, configfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    with open(configfile, 'wb') as fp:
        config_obj.write(fp)


def _check_deprecated_config_options(config_obj):
    deprecation_msgs = list()
    if 's_win_length' in config_obj or 'noise_win_length' in config_obj:
        deprecation_msgs.append(
            '> "s_win_length" and "noise_win_length" config parameters '
            'are no more\n'
            '   supported. Both are replaced by "win_length".\n'
        )
    if 'traceids' in config_obj:
        deprecation_msgs.append(
            '> "traceids" config parameter has been renamed to '
            '"traceid_mapping_file".\n'
        )
    if 'ignore_stations' in config_obj or 'use_stations' in config_obj:
        deprecation_msgs.append(
            '> "ignore_stations" and "use_stations" config parameters '
            'have been renamed to\n'
            '  "ignore_traceids" and "use_traceids", respectively.\n'
        )
    if 'dataless' in config_obj:
        deprecation_msgs.append(
            '> "dataless" config parameter has been renamed to '
            '"station_metadata".\n'
        )
    if 'clip_nmax' in config_obj:
        deprecation_msgs.append(
            '> "clip_nmax" config parameter has been renamed to '
            '"clip_max_percent".\n'
            '   Note that the new default is 5% '
            '(current value in your config file: {}%)\n'.format(
                config_obj['clip_nmax'])
        )
    if 'trace_format' in config_obj:
        deprecation_msgs.append(
            '> "trace_format" config parameter is no more supported.\n'
            '   Use "sensitivity" to manually specify how sensor sensitivity '
            'should be computed.\n'
        )
    if 'PLOT_SHOW' in config_obj:
        deprecation_msgs.append(
            '> "PLOT_SHOW" config parameter has been renamed to "plot_show".\n'
        )
    if 'PLOT_SAVE' in config_obj:
        deprecation_msgs.append(
            '> "PLOT_SAVE" config parameter has been renamed to "plot_save".\n'
        )
    if 'PLOT_SAVE_FORMAT' in config_obj:
        deprecation_msgs.append(
            '> "PLOT_SAVE_FORMAT" config parameter has been renamed to '
            '"plot_save_format".\n'
        )
    if 'vp' in config_obj:
        deprecation_msgs.append(
            '> "vp" config parameter has been renamed to "vp_source".\n'
        )
    if 'vs' in config_obj:
        deprecation_msgs.append(
            '> "vs" config parameter has been renamed to "vs_source".\n'
        )
    if 'pre_p_time' in config_obj:
        deprecation_msgs.append(
            '> "pre_p_time" config parameter has been renamed to '
            '"noise_pre_time".\n'
        )
    if 'pre_s_time' in config_obj:
        deprecation_msgs.append(
            '> "pre_s_time" config parameter has been renamed to '
            '"signal_pre_time".\n'
        )
    if 'rps_from_focal_mechanism' in config_obj:
        deprecation_msgs.append(
            '> "rps_from_focal_mechanism" config parameter has been renamed to '
            '"rp_from_focal_mechanism".\n'
        )
    if 'paz' in config_obj:
        deprecation_msgs.append(
            '> "paz" config parameter has been removed and merged with '
            '"station_metadata".\n'
        )
    if deprecation_msgs:
        sys.stderr.write(
            'Error: your config file contains deprecated parameters:\n\n')
    for msg in deprecation_msgs:
        sys.stderr.write(msg)
    if deprecation_msgs:
        sys.stderr.write(
            '\nPlease upgrade your config file manually or '
            'via the "-U" option.\n'
        )
        sys.exit(1)


def _init_plotting(plot_show):
    import matplotlib.pyplot as plt
    if not plot_show:
        plt.switch_backend('Agg')


def _check_mandatory_config_params(config_obj):
    mandatory_params = [
        'p_arrival_tolerance',
        's_arrival_tolerance',
        'noise_pre_time',
        'signal_pre_time',
        'win_length',
        'taper_halfwidth',
        'spectral_smooth_width_decades',
        'bp_freqmin_acc',
        'bp_freqmax_acc',
        'bp_freqmin_shortp',
        'bp_freqmax_shortp',
        'bp_freqmin_broadb',
        'bp_freqmax_broadb',
        'freq1_acc',
        'freq2_acc',
        'freq1_shortp',
        'freq2_shortp',
        'freq1_broadb',
        'freq2_broadb',
        'rmsmin',
        'sn_min',
        'clip_max_percent',
        'spectral_sn_min',
        'vp_source',
        'vs_source',
        'rho',
        'rpp',
        'rps',
        'geom_spread_n_exponent',
        'geom_spread_cutoff_distance',
        'f_weight',
        'weight',
        't_star_0',
        't_star_0_variability',
        'Mw_0_variability',
        'a', 'b', 'c',
        'ml_bp_freqmin',
        'ml_bp_freqmax',
        'nIQR',
        'plot_spectra_maxrows',
        'plot_traces_maxrows',
        'plot_station_text_size'
    ]
    messages = list()
    for par in mandatory_params:
        if config_obj[par] is None:
            msg = '"{}" is mandatory and cannot be None'.format(par)
            messages.append(msg)
    if messages:
        msg = '\n'.join(messages)
        sys.stderr.write(msg + '\n')
        ssp_exit(1)


def configure(options, progname, config_overrides=None):
    """
    Parse command line arguments and read config file.

    :param object options: An object containing command line options
    :param str progname: The name of the program
    :param dict config_overrides: A dictionary with parameters that override or
        extend those defined in the config file
    :return: A ``Config`` object with both command line and config options.
    """
    _check_obspy_version()
    _check_library_versions()

    configspec = _parse_configspec()
    if options.sampleconf:
        _write_sample_config(configspec, progname)
        sys.exit(0)
    if options.updateconf:
        _update_config_file(options.updateconf, configspec)
        sys.exit(0)

    config_obj = _read_config(options.config_file, configspec)

    # Apply overrides
    if config_overrides is not None:
        try:
            for key, value in config_overrides.items():
                config_obj[key] = value
        except AttributeError:
            raise ValueError('"config_override" must be a dict-like.')

    # Set to None all the 'None' strings
    for key, value in config_obj.dict().items():
        if value == 'None':
            config_obj[key] = None

    val = Validator()
    test = config_obj.validate(val)
    if isinstance(test, dict):
        for entry in test:
            if not test[entry]:
                sys.stderr.write('Invalid value for "%s": "%s"\n' %
                                 (entry, config_obj[entry]))
        sys.exit(1)
    if not test:
        sys.stderr.write('No configuration value present!\n')
        sys.exit(1)

    _check_deprecated_config_options(config_obj)
    _check_mandatory_config_params(config_obj)

    # Create a 'no_evid_' subdir into outdir.
    # The random hex string will make it sure that this name is unique
    # It will be then renamed once an evid is available
    hex = uuid.uuid4().hex
    options.outdir = os.path.join(
        options.outdir,
        'no_evid_{}'.format(hex)
    )
    _write_config(config_obj, progname, options.outdir)

    # Create a Config object
    config = Config(config_obj.dict().copy())

    # Add options to config:
    config.options = options

    # Additional config values
    config.vertical_channel_codes = ['Z']
    config.horizontal_channel_codes_1 = ['N', 'R']
    config.horizontal_channel_codes_2 = ['E', 'T']
    msc = config.mis_oriented_channels
    if msc is not None:
        config.vertical_channel_codes.append(msc[0])
        config.horizontal_channel_codes_1.append(msc[1])
        config.horizontal_channel_codes_2.append(msc[2])

    # A list of warnings to be issued when logger is set up
    config.warnings = list()

    if config.html_report:
        if not config.plot_save:
            msg = (
                'The "html_report" option is selected but "plot_save" '
                'is "False". Setting "plot_save" to "True".')
            config.warnings.append(msg)
            config.plot_save = True
        if config.plot_save_format != 'png':
            msg = (
                'The "html_report" option is selected but "plot_save_format" '
                'is not "png". Setting "plot_save_format" to "png".')
            config.warnings.append(msg)
            config.plot_save_format = 'png'
        if not config.plot_station_map:
            msg = (
                'The "html_report" option is selected but "plot_station_map" '
                'is "False". Setting "plot_station_map" to "True".')
            config.warnings.append(msg)
            config.plot_station_map = True

    if config.plot_station_map:
        try:
            _check_cartopy_version()
            _check_pyproj_version()
        except ImportError as err:
            for msg in config.warnings:
                print(msg)
            sys.stderr.write(str(err))
            sys.exit(1)

    if config.NLL_time_dir is not None or config.NLL_model_dir is not None:
        try:
            _check_nllgrid_version()
        except ImportError as err:
            sys.stderr.write(str(err))
            sys.exit(1)

    _init_plotting(config.plot_show)
    # Create a dict to store figure paths
    config.figures = defaultdict(list)
    return config


def save_config(config):
    """Save config file to output dir."""
    # Actually, it renames the file already existing.
    src = os.path.join(config.options.outdir, 'source_spec.conf')
    evid = config.hypo.evid
    dst = os.path.join(config.options.outdir, '%s.ssp.conf' % evid)
    # On Windows, dst file must not exist
    try:
        os.remove(dst)
    except Exception:
        pass
    os.rename(src, dst)


def move_outdir(config):
    """Move outdir to a new dir named from evid (and optional run_id)."""
    try:
        evid = config.hypo.evid
    except Exception:
        return
    src = config.options.outdir
    run_id = config.options.run_id
    path = os.path.normpath(config.options.outdir).split(os.sep)
    dst = os.path.join(*path[:-1], str(evid))
    if run_id:
        dst += '_' + run_id
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


def remove_old_outdir(config):
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
        if(levelno >= logging.CRITICAL):
            color = '\x1b[31;1m'  # red
        elif(levelno >= logging.ERROR):
            color = '\x1b[31;1m'  # red
        elif(levelno >= logging.WARNING):
            color = '\x1b[33;1m'  # yellow
        elif(levelno >= logging.INFO):
            # color = '\x1b[32;1m'  # green
            color = '\x1b[0m'  # no color
        elif(levelno >= logging.DEBUG):
            color = '\x1b[35;1m'  # purple
        else:
            color = '\x1b[0m'  # no color
        # Color-code the message
        args[0].msg = "{0}{1}\x1b[0m".format(color, args[0].msg)
        return fn(*args)
    return new


def setup_logging(config, basename=None, progname='source_spec'):
    """
    Set up the logging infrastructure.

    This function is typically called twice: the first time without basename
    and a second time with a basename (typically the eventid).
    When called the second time, the previous logfile is renamed using the
    given basename.
    """
    global oldlogfile
    global logger
    global PYTHON_VERSION_STR
    global OBSPY_VERSION_STR
    global NUMPY_VERSION_STR
    global SCIPY_VERSION_STR
    global MATPLOTLIB_VERSION_STR
    global CARTOPY_VERSION_STR
    # Create outdir
    if not os.path.exists(config.options.outdir):
        os.makedirs(config.options.outdir)

    if basename:
        logfile = os.path.join(config.options.outdir, '%s.ssp.log' % basename)
    else:
        datestring = datetime.now().strftime('%Y%m%d_%H%M%S')
        logfile = os.path.join(config.options.outdir,
                               '%s.ssp.log' % datestring)

    logger_root = logging.getLogger()
    if oldlogfile:
        hdlrs = logger_root.handlers[:]
        for hdlr in hdlrs:
            hdlr.flush()
            hdlr.close()
            logger_root.removeHandler(hdlr)
        # Copy old logfile to new
        shutil.copyfile(oldlogfile, logfile)
        # Remove old logfile from old and new dir
        os.remove(oldlogfile)
        oldlogfile = os.path.join(
            config.options.outdir, os.path.basename(oldlogfile))
        os.remove(oldlogfile)
        filemode = 'a'
    else:
        filemode = 'w'

    # captureWarnings is not supported in old versions of python
    try:
        logging.captureWarnings(True)
    except Exception:
        pass
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

    logger = logging.getLogger(progname)

    # Only write these debug infos for a new logfile
    if not oldlogfile:
        logger.debug('source_spec START')
        logger.debug('SourceSpec version: ' + get_versions()['version'])
        uname = platform.uname()
        uname_str = '{} {} {}'.format(uname[0], uname[2], uname[4])
        logger.debug('Platform: ' + uname_str)
        logger.debug('Python version: ' + PYTHON_VERSION_STR)
        logger.debug('ObsPy version: ' + OBSPY_VERSION_STR)
        logger.debug('NumPy version: ' + NUMPY_VERSION_STR)
        logger.debug('SciPy version: ' + SCIPY_VERSION_STR)
        logger.debug('Matplotlib version: ' + MATPLOTLIB_VERSION_STR)
        if CARTOPY_VERSION_STR is not None:
            logger.debug('Cartopy version: ' + CARTOPY_VERSION_STR)
        logger.debug('Running arguments:')
        logger.debug(' '.join(sys.argv))
    oldlogfile = logfile

    # See if there are warnings to deliver
    for _ in range(len(config.warnings)):
        msg = config.warnings.pop(0)
        logger.warning(msg)


def ssp_exit(retval=0, abort=False):
    # ssp_exit might have already been called if multiprocessing
    global ssp_exit_called
    if ssp_exit_called:
        return
    ssp_exit_called = True
    if abort:
        print('\nAborting.')
        if logger is not None:
            logger.debug('source_spec ABORTED')
    else:
        if logger is not None:
            logger.debug('source_spec END')
    logging.shutdown()
    sys.exit(retval)


def sigint_handler(sig, frame):
    ssp_exit(1, abort=True)


signal.signal(signal.SIGINT, sigint_handler)
