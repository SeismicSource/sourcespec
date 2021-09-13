# -*- coding: utf-8 -*-
"""
Setup functions for sourcespec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2021 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.utils import iteritems

import sys
import os
import shutil
import logging
import signal
import numpy as np
from argparse import ArgumentParser, RawTextHelpFormatter
from datetime import datetime
try:
    from itertools import izip_longest as zip_longest
except ImportError:
    from itertools import zip_longest
from sourcespec.configobj import ConfigObj
from sourcespec.configobj.validate import Validator
from sourcespec.config import Config
from sourcespec.AsyncPlotter import AsyncPlotter
from sourcespec._version import get_versions

# define ipshell(), if possible
if sys.stdout.isatty():
    try:
        import IPython
        IP_VERSION = tuple(map(int, IPython.__version__.split('.')[:3]))
        len(IP_VERSION) > 2 or IP_VERSION.append(0)
        if IP_VERSION < (0, 11, 0):
            from IPython.Shell import IPShellEmbed
            ipshell = IPShellEmbed()
        if (0, 11, 0) <= IP_VERSION < (1, 0, 0):
            from IPython.frontend.terminal.embed import InteractiveShellEmbed
            ipshell = InteractiveShellEmbed()
        elif IP_VERSION >= (1, 0, 0):
            from IPython.terminal.embed import InteractiveShellEmbed
            ipshell = InteractiveShellEmbed()
    except ImportError:
        ipshell = None
else:
    ipshell = None

# global variables
OS = os.name
DEBUG = False
oldlogfile = None
plotter = None
logger = None
ssp_exit_called = False
OBSPY_VERSION = None
OBSPY_VERSION_STR = None
NUMPY_VERSION_STR = None
SCIPY_VERSION_STR = None
MATPLOTLIB_VERSION_STR = None
PYTHON_VERSION_STR = None


def _check_obspy_version():
    global OBSPY_VERSION, OBSPY_VERSION_STR
    # check ObsPy version
    import obspy
    MIN_OBSPY_VERSION = (1, 0, 0)
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


def dprint(string):
    """Print debug information."""
    if DEBUG:
        sys.stderr.write(string)
        sys.stderr.write('\n')


def _parse_values(value_str):
    if value_str[0] == 'i':
        value_str = value_str[1:]
        val_min, val_max, val_step =\
            tuple(map(float, value_str.rstrip(',').split(',')))
        # we add a small number to val_max to make sure
        # it is included by np.arange()
        output = tuple(np.arange(val_min, val_max + 1e-9, val_step))
    else:
        try:
            output = tuple(map(float, value_str.rstrip(',').split(',')))
        except ValueError:
            sys.stderr.write('ERROR: Invalid value: %s\n' % value_str)
            sys.exit(1)
    return output


def _get_description(progname):
    if progname == 'source_spec':
        nargs = '+'
        description = 'Estimation of seismic source parameters from '
        description += 'inversion of S-wave spectra.'
        epilog = ''
    elif progname == 'source_model':
        nargs = '*'
        description = 'Direct modeling of S-wave spectra.'
        epilog = 'Several values of moment magnitude, seismic moment, t-star\n'
        epilog += 'and alpha can be specified using a comma-separated list,'
        epilog += 'eg.:\n'
        epilog += '  --mag=3,3.5,4,4.5\n\n'
        epilog += 'A value interval can be specified by prepending "i" and\n'
        epilog += 'indicating min, max and step, eg.:\n'
        epilog += '  --fc=i1.0,5.0,0.5\n\n'
        epilog += 'When specifing several values, by default a simple\n'
        epilog += 'correspondance is established, e.g.:\n'
        epilog += '  --mag=3,3.5 --fc=1.0,2.0,3.0\n'
        epilog += 'will generate the couples:\n'
        epilog += '  (3, 1.0), (3.5, 2.0), (3.5, 3.0)\n'
        epilog += '(note that the magnitude value "3.5" is repeated twice).\n'
        epilog += 'Use "-C" to generate all the possible combinations.'
    else:
        sys.stderr.write('Wrong program name: %s\n' % progname)
        sys.exit(1)
    return description, epilog, nargs


def _init_parser(description, epilog, nargs):
    parser = ArgumentParser(description=description,
                            epilog=epilog,
                            formatter_class=RawTextHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-S', '--sampleconf', dest='sampleconf',
                       action='store_true', default=False,
                       help='write sample configuration to file and exit')
    group.add_argument('-U', '--updateconf', dest='updateconf',
                       action='store', default=None,
                       help='update an existing config file from a previous '
                       'version', metavar='FILE')
    parser.add_argument('-c', '--configfile', dest='config_file',
                        action='store', default='source_spec.conf',
                        help='load configuration from FILE '
                        '(default: source_spec.conf)', metavar='FILE')
    group.add_argument('-t', '--trace_path', nargs=nargs,
                       help='path to trace file(s) or trace dir')
    parser.add_argument('-H', '--hypocenter', dest='hypo_file',
                        action='store', default=None,
                        help='get hypocenter information from FILE',
                        metavar='FILE')
    parser.add_argument('-p', '--pickfile', dest='pick_file',
                        action='store', default=None,
                        help='get picks from FILE', metavar='FILE')
    parser.add_argument('-q', '--qmlfile', dest='qml_file',
                        action='store', default=None,
                        help='get picks and hypocenter information from '
                        'QuakeML FILE',
                        metavar='FILE')
    parser.add_argument('-n', '--evname', dest='evname',
                        action='store', default=None,
                        help='event name (used for plots and output files) ',
                        metavar='NAME')
    parser.add_argument('-e', '--evid', dest='evid',
                        action='store', default=None,
                        help='get evid from catalog', metavar='EVID')
    parser.add_argument('-s', '--station', dest='station',
                        action='store', default=None,
                        help='only use this station', metavar='STATION')
    parser.add_argument('-N', '--no-response', dest='no_response',
                        action='store_true', default=False,
                        help='do not remove instrument response')
    return parser


def _update_parser(parser, progname):
    if progname == 'source_spec':
        parser.add_argument('-o', '--outdir', dest='outdir',
                            action='store', default='sspec_out',
                            help='save output to OUTDIR (default: sspec_out)',
                            metavar='OUTDIR')
        parser.add_argument('-C', '--correction', dest='correction',
                            action='store_true', default=False,
                            help='apply station correction to the "H" '
                            'component of the spectra')
    elif progname == 'source_model':
        parser.add_argument('-f', '--fmin', dest='fmin', action='store',
                            type=float, default='0.01',
                            help='minimum frequency (Hz, default 0.01)',
                            metavar='FMIN')
        parser.add_argument('-F', '--fmax', dest='fmax', action='store',
                            type=float, default='50.0',
                            help='maximum frequency (Hz, default 50.0)',
                            metavar='FMAX')
        parser.add_argument('-k', '--fc', dest='fc', action='store',
                            default='10.0',
                            help='(list of) corner frequency '
                            '(Hz, default 10.0)', metavar='FC')
        parser.add_argument('-m', '--mag', dest='mag', action='store',
                            default='2.0',
                            help='(list of) moment magnitude (default 2.0)',
                            metavar='Mw')
        parser.add_argument('-M', '--moment', dest='Mo', action='store',
                            default='NaN',
                            help='(list of) seismic moment '
                            '(N.m, default undefined)', metavar='Mo')
        parser.add_argument('-*', '--tstar', dest='t_star', action='store',
                            default='0.0',
                            help='(list of) t-star (attenuation, default 0.0)',
                            metavar='T-STAR')
        parser.add_argument('-a', '--alpha', dest='alpha', action='store',
                            default='1.0',
                            help='(list of) alpha (exponent for frequency '
                            'dependence\nof attenuation, default 1.0)',
                            metavar='1.0')
        parser.add_argument('-C', '--combine', dest='combine',
                            action='store_true', default=False,
                            help='generate all the combinations of fc, mag, '
                            'Mo, tstar, alpha')
        parser.add_argument('-P', '--plot', dest='plot', action='store_true',
                            default=False, help='plot results')
    parser.add_argument('-v', '--version', action='version',
                        version=get_versions()['version'])


def _parse_args(progname):
    """Parse command line arguments."""
    description, epilog, nargs = _get_description(progname)
    parser = _init_parser(description, epilog, nargs)
    _update_parser(parser, progname)

    options = parser.parse_args()

    if options.sampleconf:
        configspec = _parse_configspec()
        _write_sample_config(configspec, 'source_spec')
        sys.exit(0)

    if options.updateconf:
        _update_config_file(options.updateconf)
        sys.exit(0)

    if progname == 'source_model':
        options.mag = _parse_values(options.mag)
        options.Mo = _parse_values(options.Mo)
        options.alpha = _parse_values(options.alpha)
        options.fc = _parse_values(options.fc)
        options.t_star = _parse_values(options.t_star)

        if options.combine:
            oplist = [(fc, mag, Mo, t_star, alpha)
                      for fc in options.fc
                      for mag in options.mag
                      for Mo in options.Mo
                      for t_star in options.t_star
                      for alpha in options.alpha]
            oplist = list(map(list, zip(*oplist)))
        else:
            # Add trailing "None" to shorter lists and zip:
            oplist = [options.fc, options.mag, options.Mo,
                      options.t_star, options.alpha]
            oplist = list(zip_longest(*oplist))
            # Unzip and convert tuple to lists:
            oplist = list(map(list, zip(*oplist)))
            for l in oplist:
                for n, x in enumerate(l):
                    if x is None:
                        l[n] = l[n-1]

        options.fc, options.mag, options.Mo, options.t_star, options.alpha = \
            oplist
        # Add unused options (required by source_spec):
        options.correction = False
        options.outdir = '.'

    return options


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
    configfile = progname + '.conf'
    write_file = True
    if os.path.exists(configfile):
        # Workaround for python2/3 compatibility
        try:
            raw_input = input
        except NameError:
            pass
        ans = raw_input('%s already exists. '
                        'Do you want to overwrite it? [y/N] ' % configfile)
        if ans in ['y', 'Y']:
            write_file = True
        else:
            write_file = False
    if write_file:
        with open(configfile, 'wb') as fp:
            c.write(fp)
        print('Sample config file written to: ' + configfile)


def _update_config_file(config_file):
    configspec = _parse_configspec()
    config_obj = _read_config(config_file, configspec)
    val = Validator()
    config_obj.validate(val)
    mod_time = datetime.fromtimestamp(os.path.getmtime(config_file))
    mod_time_str = mod_time.strftime('%Y%m%d_%H%M%S')
    config_file_old = '{}.{}'.format(config_file, mod_time_str)
    # Workaround for python2/3 compatibility
    try:
        raw_input = input
    except NameError:
        pass
    ans = raw_input(
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
    for k, v in config_obj.items():
        if k not in config_new:
            continue
        config_new[k] = v
    # Migrate 's_win_length' to 'win_length'
    if 's_win_length' in config_obj:
        config_new['win_length'] = config_obj['s_win_length']
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


def configure(progname='source_spec'):
    """
    Parse command line arguments and read config file.

    Returns a ``Config`` object.
    """
    _check_obspy_version()
    _check_library_versions()

    global DEBUG
    options = _parse_args(progname)

    configspec = _parse_configspec()
    config_obj = _read_config(options.config_file, configspec)

    # Set to None all the 'None' strings
    for key, value in iteritems(config_obj.dict()):
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

    if 's_win_length' in config_obj or 'noise_win_length' in config_obj:
        sys.stderr.write(
            'Error: "s_win_length" and "noise_win_length" config parameters '
            'are no more\nsupported. Both are replaced by "win_length".\n\n'
            'Consider upgrading your config file via the "-U" option.\n'
        )
        sys.exit(1)

    _write_config(config_obj, progname, options.outdir)

    # Create a Config object
    config = Config(config_obj.dict().copy())
    DEBUG = config.DEBUG

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

    if config.plot_station_map:
        _exit = False
        cartopy_min_ver = (0, 18, 0)
        try:
            cartopy_ver = None
            import cartopy  #NOQA
            cartopy_ver = tuple(map(int, cartopy.__version__.split('.')[:3]))
            if cartopy_ver < cartopy_min_ver:
                raise ImportError
        except ImportError:
            sys.stderr.write(
                'Please install cartopy >= {}.{}.{} to plot maps.\n'.format(
                    *cartopy_min_ver)
            )
            if cartopy_ver:
                sys.stderr.write(
                    'Installed cartopy version: {}.{}.{}.\n'.format(
                        *cartopy_ver)
                )
            _exit = True
        try:
            import pyproj  #NOQA
        except ImportError:
            sys.stderr.write('Please install pyproj to plot maps.\n')
            _exit = True
        if _exit:
            sys.exit(1)

    return config


def save_config(config):
    """Save config file to output dir."""
    # Actually, it renames the file already existing.
    src = os.path.join(config.options.outdir, 'source_spec.conf')
    evid = config.hypo.evid
    dst = os.path.join(config.options.outdir, '%s.ssp.conf' % evid)
    os.rename(src, dst)


def setup_logging(config, basename=None, progname='source_spec'):
    """
    Setup the logging infrastructure.

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
        if OS == 'nt':
            # Windows has problems with renaming
            shutil.copyfile(oldlogfile, logfile)
        else:
            os.rename(oldlogfile, logfile)
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
    logger_root.addHandler(console)

    logger = logging.getLogger(progname)

    # Only write these debug infos for a new logfile
    if not oldlogfile:
        logger.debug('source_spec START')
        logger.debug('SourceSpec version: ' + get_versions()['version'])
        logger.debug('Python version: ' + PYTHON_VERSION_STR)
        logger.debug('ObsPy version: ' + OBSPY_VERSION_STR)
        logger.debug('NumPy version: ' + NUMPY_VERSION_STR)
        logger.debug('SciPy version: ' + SCIPY_VERSION_STR)
        logger.debug('Matplotlib version: ' + MATPLOTLIB_VERSION_STR)
        logger.debug('Running arguments:')
        logger.debug(' '.join(sys.argv))
    oldlogfile = logfile


def init_plotting(config):
    if config.html_report:
        if not config.PLOT_SAVE:
            logger.warning(
                'The html_report option is selected but PLOT_SAVE is False. '
                'Setting PLOT_SAVE to True.')
            config.PLOT_SAVE = True
        if config.PLOT_SAVE_FORMAT != 'png':
            logger.warning(
                'The html_report option is selected but PLOT_SAVE_FORMAT is '
                'not png. Setting PLOT_SAVE_FORMAT to png.')
            config.PLOT_SAVE_FORMAT = 'png'
        if not config.plot_station_map:
            logger.warning(
                'The html_report option is selected but plot_station_map '
                'is False. Setting plot_station_map to True.')
            config.plot_station_map = True

    import matplotlib.pyplot as plt
    if not config.PLOT_SHOW:
        plt.switch_backend('Agg')
    global plotter
    if OS == 'nt':
        # AsyncPlotter() doesn't work on Windows. TODO: fix this
        plotter = None
    else:
        plotter = AsyncPlotter()
    return plotter


def ssp_exit(retval=0, abort=False):
    # ssp_exit might have already been called if multiprocessing
    global ssp_exit_called
    if ssp_exit_called:
        return
    ssp_exit_called = True
    if abort:
        print('\nAborting.')
        logger.debug('source_spec ABORTED')
    else:
        logger.debug('source_spec END')
    logging.shutdown()
    global plotter
    if plotter is not None:
        if abort:
            plotter.terminate()
        else:
            plotter.join()
    sys.exit(retval)


def sigint_handler(sig, frame):
    ssp_exit(1, abort=True)


signal.signal(signal.SIGINT, sigint_handler)
