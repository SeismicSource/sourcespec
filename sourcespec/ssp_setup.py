# -*- coding: utf-8 -*-
# ssp_setup.py
#
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013-2014 Claudio Satriano <satriano@ipgp.fr>,
#               Emanuela Matrullo <matrullo@geologie.ens.fr>,
#               Agnes Chounet <chounet@ipgp.fr>
# (c) 2015-2016 Claudio Satriano <satriano@ipgp.fr>
"""Setup functions for source_spec."""
import sys
import os
import shutil
import logging
import signal
import numpy as np
from configobj import ConfigObj
from validate import Validator
from config import Config
from argparse import ArgumentParser, RawTextHelpFormatter
from version import get_git_version
from datetime import datetime
from multiprocessing import Pool, cpu_count

# check ObsPy version
import obspy
MIN_OBSPY_VERSION = (1, 0, 0)
OBSPY_VERSION_STR = obspy.__version__
OBSPY_VERSION = tuple(map(int, OBSPY_VERSION_STR.split('.')[:3]))
if OBSPY_VERSION < MIN_OBSPY_VERSION:
    msg = 'ERROR: ObsPy >= %s.%s.%s is required.' % MIN_OBSPY_VERSION
    msg += ' You have version: %s\n' % OBSPY_VERSION_STR
    sys.stderr.write(msg)
    sys.exit(1)

# define ipshell(), if possible
if sys.stdout.isatty():
    try:
        import IPython
        IP_VERSION = map(int, IPython.__version__.split('.')[:3])
        len(IP_VERSION) > 2 or IP_VERSION.append(0)
        if IP_VERSION < [0, 11, 0]:
            from IPython.Shell import IPShellEmbed
            ipshell = IPShellEmbed()
        if [0, 11, 0] <= IP_VERSION < [1, 0, 0]:
            from IPython.frontend.terminal.embed import InteractiveShellEmbed
            ipshell = InteractiveShellEmbed()
        elif IP_VERSION >= [1, 0, 0]:
            from IPython.terminal.embed import InteractiveShellEmbed
            ipshell = InteractiveShellEmbed()
    except ImportError:
        ipshell = None
else:
    ipshell = None

# global variables
OS = os.name
DEBUG = False
matplotlib_unloaded = False
oldlogfile = None
plot_pool = None
ssp_exit_called = False


def dprint(string):
    """Print debug information."""
    if DEBUG:
        sys.stderr.write(string)
        sys.stderr.write('\n')


def unload_matplotlib():
    global matplotlib_unloaded
    if matplotlib_unloaded:
        return
    # Source:
    #   http://stackoverflow.com/q/3285193
    modules = []
    for module in sys.modules:
        if module.startswith('matplotlib'):
            modules.append(module)
    for module in modules:
        sys.modules.pop(module)
    matplotlib_unloaded = True


def _parse_values(value_str):
    if value_str[0] == 'i':
        value_str = value_str[1:]
        val_min, val_max, val_step =\
            map(float, value_str.rstrip(',').split(','))
        # we add a small number to val_max to make sure
        # it is included by np.arange()
        output = tuple(np.arange(val_min, val_max + 1e-9, val_step))
    else:
        try:
            output = map(float, value_str.rstrip(',').split(','))
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
    parser.add_argument('-c', '--configfile', dest='config_file',
                        action='store', default='config.py',
                        help='load configuration from FILE '
                        '(default: config.py)', metavar='FILE')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-t', '--trace_path', nargs=nargs,
                       help='path to trace file(s) or trace dir')
    group.add_argument('-S', '--sampleconf', dest='sampleconf',
                       action='store_true', default=False,
                       help='write sample configuration to file and exit')
    parser.add_argument('-H', '--hypocenter', dest='hypo_file',
                        action='store', default=None,
                        help='get hypocenter information from FILE',
                        metavar='FILE')
    parser.add_argument('-p', '--pickfile', dest='pick_file',
                        action='store', default=None,
                        help='get picks from FILE', metavar='FILE')
    parser.add_argument('-e', '--evid', dest='evid',
                        action='store', default=None,
                        help='get evid from catalog', metavar='EVID')
    parser.add_argument('-s', '--station', dest='station',
                        action='store', default=None,
                        help='only use this station', metavar='STATION')
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
                        version=get_git_version())


def _parse_args(progname):
    """"Parse command line arguments."""
    description, epilog, nargs = _get_description(progname)
    parser = _init_parser(description, epilog, nargs)
    _update_parser(parser, progname)

    options = parser.parse_args()

    if options.sampleconf:
        configspec = _parse_configspec()
        _write_sample_config(configspec, 'source_spec')
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
            oplist = map(list, zip(*oplist))
        else:
            # Add trailing "None" to shorter lists and zip:
            oplist = [options.fc, options.mag, options.Mo,
                      options.t_star, options.alpha]
            oplist = map(None, *oplist)
            # Unzip and convert tuple to lists:
            oplist = map(list, zip(*oplist))
            for l in oplist:
                for n, x in enumerate(l):
                    if x is None:
                        l[n] = l[n-1]

        options.fc, options.mag, options.Mo, options.t_star, options.alpha = \
            oplist
        # Add unused options (required by source_spec):
        options.correction = False

    return options


def _parse_configspec():
    try:
        configspec_file = os.path.join(os.path.dirname(__file__),
                                       'configspec.conf')
        configspec = ConfigObj(configspec_file, interpolation=False,
                               list_values=False, _inspec=True,
                               file_error=True)
    except IOError, message:
        sys.stderr.write('%s\n' % message)
        sys.exit(1)
    except Exception, message:
        sys.stderr.write('Unable to read "%s": %s\n' %
                         (configspec_file, message))
        sys.exit(1)
    return configspec


def _write_sample_config(configspec, progname):
    c = ConfigObj(configspec=configspec)
    val = Validator()
    c.validate(val)
    c.defaults = []
    c.initial_comment = configspec.initial_comment
    c.comments = configspec.comments
    configfile = progname + '.conf'
    with open(configfile, 'w') as fp:
        c.write(fp)
    print 'Sample config file written to: ' + configfile


def configure(progname='source_spec'):
    """
    Parse command line arguments and read config file.

    Returns a ``Config`` object.
    """
    global DEBUG
    options = _parse_args(progname)

    configspec = _parse_configspec()
    try:
        config_file = options.config_file
        config_obj = ConfigObj(config_file, configspec=configspec,
                               file_error=True)
    except IOError, message:
        sys.stderr.write('%s\n' % message)
        sys.exit(1)
    except Exception, message:
        sys.stderr.write('Unable to read "%s": %s\n' % (config_file, message))
        sys.exit(1)

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

    # Create a Config object
    config = Config(config_obj.dict().copy())
    DEBUG = config.DEBUG

    #add options to config:
    config.options = options

    return config


def setup_logging(config, basename=None):
    """Setup the logging infrastructure."""
    global oldlogfile
    # Create outdir
    if not os.path.exists(config.options.outdir):
        os.makedirs(config.options.outdir)

    if basename:
        logfile = os.path.join(config.options.outdir, '%s.ssp.log' % basename)
    else:
        datestring = datetime.now().strftime('%Y%m%d_%H%M%S')
        logfile = os.path.join(config.options.outdir,
                               '%s.ssp.log' % datestring)

    log = logging.getLogger()
    if oldlogfile:
        hdlrs = log.handlers[:]
        for hdlr in hdlrs:
            hdlr.flush()
            hdlr.close()
            log.removeHandler(hdlr)
        if OS == 'nt':
            # Windows has problems with renaming
            shutil.copyfile(oldlogfile, logfile)
        else:
            os.rename(oldlogfile, logfile)
        filemode = 'a'
    else:
        filemode = 'w'
    oldlogfile = logfile

    #logging.basicConfig(
    #        level=logging.DEBUG,
    #        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
    #        filename=logfile,
    #        filemode='w'
    #        )

    # captureWarnings is not supported in old versions of python
    try:
        logging.captureWarnings(True)
    except:
        pass
    log.setLevel(logging.DEBUG)
    filehand = logging.FileHandler(filename=logfile, mode=filemode)
    filehand.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s %(name)-12s '
                                  '%(levelname)-8s %(message)s')
    filehand.setFormatter(formatter)
    log.addHandler(filehand)

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    log.addHandler(console)

    if not basename:
        logging.debug('source_spec START')
        logging.debug('Version: ' + get_git_version())
        # write running arguments
        logging.debug('Running arguments:')
        logging.debug(' '.join(sys.argv))


def init_plotting():
    global plot_pool
    if plot_pool is None:
        processes = cpu_count() - 1
        processes = processes if processes > 0 else 1
        plot_pool = Pool(processes)
    return plot_pool


def ssp_exit(retval=0, message=''):
    # ssp_exit might have already been called if multiprocessing
    global ssp_exit_called
    if ssp_exit_called:
        return
    ssp_exit_called = True
    if message:
        logging.debug(message)
    else:
        logging.debug('source_spec END')
    logging.shutdown()
    global plot_pool
    if plot_pool is not None:
        plot_pool.close()
        plot_pool.join()
    sys.exit(retval)


def sigint_handler(sig, frame):
    print ''
    print 'Aborting.'
    ssp_exit(1, 'source_spec ABORTED')
signal.signal(signal.SIGINT, sigint_handler)
