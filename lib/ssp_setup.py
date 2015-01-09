# -*- coding: utf-8 -*-
# ssp_setup.py
#
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013-2014 Claudio Satriano <satriano@ipgp.fr>,
#               Emanuela Matrullo <matrullo@geologie.ens.fr>,
#               Agnes Chounet <chounet@ipgp.fr>
# (c) 2015 Claudio Satriano <satriano@ipgp.fr>
'''
Setup functions for source_spec
'''
import sys
import os
import logging
import numpy as np
from configobj import ConfigObj
from validate import Validator
from config import Config
from argparse import ArgumentParser
from datetime import datetime


if sys.stdout.isatty():
    import IPython
    ip_version = map(int, IPython.__version__.split('.'))
    if ip_version[0] == 0:
        if ip_version[1] >= 11:
        # ipython >= 0.11
            try:
                from IPython.frontend.terminal.embed import InteractiveShellEmbed
                ipshell = InteractiveShellEmbed()
            except ImportError:
                ipshell = None
        else:
        # ipython < 0.11
            try:
                from IPython.Shell import IPShellEmbed
                ipshell = IPShellEmbed()
            except ImportError:
                ipshell = None
    elif ip_version[0] > 0:
    # ipython >= 1.0.0
        try:
            from IPython.terminal.embed import InteractiveShellEmbed
            ipshell = InteractiveShellEmbed()
        except ImportError:
            ipshell = None
else:
    ipshell = None


DEBUG = False
def dprint(string):
    '''
    Print debug information.
    '''
    if DEBUG:
        sys.stderr.write(string)
        sys.stderr.write('\n')


matplotlib_unloaded = False
def unload_matplotlib():
    global matplotlib_unloaded
    if matplotlib_unloaded:
        return
    # Source:
    #   http://stackoverflow.com/questions/3285193/how-to-switch-backends-in-matplotlib-python
    modules = []
    for module in sys.modules:
        if module.startswith('matplotlib'):
            modules.append(module)
    for module in modules:
        sys.modules.pop(module)
    matplotlib_unloaded = True


def __parse_args(progname):
    '''
    Parse command line arguments
    '''
    if progname == 'source_spec':
        nargs='+'
    elif progname == 'source_model':
        nargs='*'
    else:
        sys.stderr.write('Wrong program name: %s\n' % progname)
        sys.exit(1)

    parser = ArgumentParser()
    parser.add_argument('-c', '--configfile', dest='config_file', action='store', default='config.py',
            help='load configuration from FILE (default: config.py)', metavar='FILE')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-t', '--trace_path', nargs=nargs, help='path to trace file(s) or trace dir')
    group.add_argument('-S', '--sampleconf', dest='sampleconf', action='store_true', default=False,
            help='write sample configuration to file and exit')
    parser.add_argument('-H', '--hypocenter', dest='hypo_file', action='store', default=None,
            help='get hypocenter information from FILE', metavar='FILE')
    parser.add_argument('-p', '--pickfile', dest='pick_file', action='store', default=None,
            help='get picks from FILE', metavar='FILE')
    parser.add_argument('-e', '--evid', dest='evid', action='store', default=None,
            help='get evid from catalog', metavar='EVID')
    parser.add_argument('-s', '--station', dest='station', action='store', default=None,
            help='only use this station', metavar='STATION')
    if progname == 'source_spec':
        parser.add_argument('-o', '--outdir', dest='outdir', action='store', default='sspec_out',
                help='save output to OUTDIR (default: sspec_out)', metavar='OUTDIR')
        parser.add_argument('-C', '--correction', dest='correction', action='store_true', default=False,
                help='apply station correction to the "H" component of the spectra')
    elif progname == 'source_model':
        parser.add_argument('-f', '--fmin', dest='fmin', action='store', default='0.01',
                help='minimum frequency (Hz, default 0.01)', metavar='FMIN')
        parser.add_argument('-F', '--fmax', dest='fmax', action='store', default='50.0',
                help='maximum frequency (Hz, default 50.0)', metavar='FMAX')
        parser.add_argument('-k', '--fc', dest='fc', action='store', default='10.0',
                help='corner frequency (Hz, default 10.0)', metavar='FC')
        parser.add_argument('-m', '--mag', dest='mag', action='store', default='2.0',
                help='moment magnitude (default 2.0)', metavar='Mw')
        parser.add_argument('-M', '--moment', dest='Mo', action='store', default='NaN',
                help='seismic moment (N.m, default undefined)', metavar='Mo')
        parser.add_argument('-*', '--tstar', dest='t_star', action='store', default='0.0',
                help='t-star (attenuation, default 0.0)', metavar='T-STAR')
        parser.add_argument('-a', '--alpha', dest='alpha', action='store', default='1.0',
                help='alpha (exponent for frequency dependence of attenuation, default 1.0)', metavar='1.0')
        parser.add_argument('-C', '--combine', dest='combine', action='store_true', default=False,
                help='generate all the combinations of fc, mag, Mo, tstar')
        parser.add_argument('-P', '--plot', dest='plot', action='store_true', default=False,
                help='plot results')

    options = parser.parse_args()

    if options.sampleconf:
        configspec = __parse_configspec()
        __write_sample_config(configspec, 'source_spec')
        sys.exit(0)

    if progname == 'source_model':
        options.fmin = float(options.fmin)
        options.fmax = float(options.fmax)

        options.mag = map(float, options.mag.rstrip(',').split(','))
        options.Mo = map(float, options.Mo.rstrip(',').split(','))
        options.alpha = map(float, options.alpha.rstrip(',').split(','))
        if options.fc[0] == 'i' or options.t_star[0] == 'i':
            if options.fc[0] == 'i':
                options.fc = options.fc[1:]
                fc_min, fc_max, fc_step = map(float, options.fc.rstrip(',').split(','))
                options.fc = tuple(np.arange(fc_min, fc_max+fc_step, fc_step))
            else:
                options.fc = map(float, options.fc.rstrip(',').split(','))

            if options.t_star[0] == 'i':
                options.t_star = options.t_star[1:]
                t_star_min, t_star_max, t_star_step = map(float, options.t_star.rstrip(',').split(','))
                options.t_star = tuple(np.arange(t_star_min, t_star_max+t_star_step, t_star_step))
            else:
                options.t_star = map(float, options.t_star.rstrip(',').split(','))

        else:
            options.fc = map(float, options.fc.rstrip(',').split(','))
            options.t_star = map(float, options.t_star.rstrip(',').split(','))

        if options.combine:
            oplist = [(fc, mag, Mo, t_star, alpha)
                    for fc in options.fc
                    for mag in options.mag
                    for Mo in options.Mo
                    for t_star in options.t_star
                    for alpha in options.alpha
                    ]
            oplist = map(list, zip(*oplist))
        else:
            # Add trailing "None" to shorter lists and zip:
            oplist = [options.fc, options.mag, options.Mo, options.t_star, options.alpha]
            oplist = map(None, *oplist)
            # Unzip and convert tuple to lists:
            oplist = map(list, zip(*oplist))
            for l in oplist:
                for n, x in enumerate(l):
                    if x is None:
                        l[n] = l[n-1]

        options.fc, options.mag, options.Mo, options.t_star, options.alpha = oplist
        # Add unused options (required by source_spec):
        options.pick_file = None
        options.correction = False

    return options


def __parse_configspec():
    try:
        configspec_file = os.path.join(os.path.dirname(__file__), 'configspec.conf')
        configspec = ConfigObj(configspec_file, interpolation=False,
                        list_values=False, _inspec=True, file_error=True)
    except IOError, message:
        sys.stderr.write('%s\n' % message)
        sys.exit(1)
    except Exception, message:
        sys.stderr.write('Unable to read "%s": %s\n' % (configspec_file, message))
        sys.exit(1)
    return configspec


def __write_sample_config(configspec, progname):
    c = ConfigObj(configspec=configspec)
    val = Validator()
    c.validate(val)
    c.defaults=[]
    c.initial_comment = configspec.initial_comment
    c.comments = configspec.comments
    configfile = progname + '.conf'
    with open(configfile, 'w') as fp:
        c.write(fp)
    print 'Sample config file written to: ' + configfile


def configure(progname='source_spec'):
    '''
    Parse command line arguments and read config file.
    Returns a ``Config`` object.
    '''
    global DEBUG
    options = __parse_args(progname)

    configspec = __parse_configspec()
    try:
        config_file = options.config_file
        config_obj = ConfigObj(config_file, configspec=configspec, file_error=True)
    except IOError, message:
        sys.stderr.write('%s\n' % message)
        sys.exit(1)
    except Exception, message:
        sys.stderr.write('Unable to read "%s": %s\n' % (config_file, message))
        sys.exit(1)

    val = Validator()
    test = config_obj.validate(val)
    if test is True:
        # no problem
        pass
    elif test is False:
        sys.stderr.write('No configuration value present!\n')
        sys.exit(1)
    else:
        for entry in test:
            if test[entry] == False:
                sys.stderr.write('Invalid value for "%s": "%s"\n'
                        % (entry, config_obj[entry]))
        sys.exit(1)

    # Create a Config object
    config = Config(config_obj.dict().copy())
    DEBUG  = config.DEBUG

    #add options to config:
    config.options = options

    return config


oldlogfile = None
def setup_logging(config, basename=None):
    '''
    Setup the logging infrastructure.
    '''
    global oldlogfile
    # Create outdir
    if not os.path.exists(config.options.outdir):
        os.makedirs(config.options.outdir)

    if basename:
        logfile = os.path.join(config.options.outdir, '%s.ssp.log' % basename)
    else:
        datestring = datetime.now().strftime('%Y%m%d_%H%M%S')
        logfile = os.path.join(config.options.outdir, '%s.ssp.log' % datestring)

    log=logging.getLogger()
    if oldlogfile:
        hdlrs = log.handlers[:]
        for hdlr in hdlrs:
            log.removeHandler(hdlr)
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
    try: logging.captureWarnings(True)
    except: pass
    log.setLevel(logging.DEBUG)
    filehand = logging.FileHandler(filename=logfile, mode=filemode)
    filehand.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    filehand.setFormatter(formatter)
    log.addHandler(filehand)

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    log.addHandler(console)

    if not basename:
        logging.debug('source_spec START')
        # write running arguments
        logging.debug('Running arguments:')
        logging.debug(' '.join(sys.argv))


def ssp_exit(retval=0):
    logging.debug('source_spec END')
    logging.shutdown()
    sys.exit(retval)
