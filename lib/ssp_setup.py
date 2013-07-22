# -*- coding: utf-8 -*-
# ssp_setup.py
#
# Setup functions for source_spec
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013 Claudio Satriano <satriano@ipgp.fr>,
#          Emanuela Matrullo <matrullo@geologie.ens.fr>,
#          Agnes Chounet <chounet@ipgp.fr>
import sys
import os
import logging
import numpy as np
from configobj import ConfigObj
from validate import Validator
from config import Config
from optparse import OptionParser
from datetime import datetime

if sys.stdout.isatty():
    try:
        # ipython >= 0.11
        from IPython.frontend.terminal.embed import InteractiveShellEmbed
        ipshell = InteractiveShellEmbed()
    except ImportError:
        try:
            # ipython < 0.11
            from IPython.Shell import IPShellEmbed
            ipshell = IPShellEmbed()
        except ImportError:
            # Give up!
            ipshell = None
else:
    ipshell = None

DEBUG=False
def dprint(string):
    if DEBUG:
        sys.stderr.write(string)
        sys.stderr.write('\n')

def __parse_args_source_spec():
    usage = 'usage: %prog [options] trace_file(s) | trace_dir'

    parser = OptionParser(usage=usage);
    parser.add_option('-c', '--configfile', dest='config_file', action='store', default='config.py',
            help='Load configuration from FILE (default: config.py)', metavar='DIR | FILE')
    parser.add_option('-H', '--hypocenter', dest='hypo_file', action='store', default=None,
            help='Get hypocenter information from FILE', metavar='FILE')
    parser.add_option('-p', '--pickfile', dest='pick_file', action='store', default=None,
            help='Get picks from FILE', metavar='FILE')
    parser.add_option('-e', '--evid', dest='evid', action='store', default=None,
            help='Get evid from catalog', metavar='EVID')
    parser.add_option('-o', '--outdir', dest='outdir', action='store', default='sspec_out',
            help='Save output to OUTDIR (default: sspec_out)', metavar='OUTDIR')
    parser.add_option('-s', '--station', dest='station', action='store', default=None,
            help='Only use this station', metavar='STATION')
    parser.add_option('-C', '--correction', dest='correction', action='store_true', default=False,
            help='Apply station correction to the "H" component of the spectra')
    parser.add_option('-S', '--sampleconf', dest='sampleconf', action='store_true', default=False,
            help='Write sample configuration to file and exit')
    (options, args) = parser.parse_args()

    if options.sampleconf:
        configspec = __parse_configspec()
        __write_sample_config(configspec, 'source_spec')
        sys.exit(0)

    if len(args) < 1:
        parser.print_usage(file=sys.stderr)
        sys.stderr.write("\tUse '-h' for help\n\n")
        sys.exit(1)

    return options, args

def __parse_args_source_model():
    usage = 'usage: %prog [options] trace_file(s) | trace_dir'

    parser = OptionParser(usage=usage);
    parser.add_option('-c', '--configfile', dest='config_file', action='store', default='config.py',
            help='Load configuration from FILE (default: config.py)', metavar='DIR | FILE')
    parser.add_option('-H', '--hypocenter', dest='hypo_file', action='store', default=None,
            help='Get hypocenter information from FILE', metavar='FILE')
    parser.add_option('-e', '--evid', dest='evid', action='store', default=None,
            help='Get evid from catalog', metavar='EVID')
    parser.add_option('-s', '--station', dest='station', action='store', default=None,
            help='Only use this station', metavar='STATION')
    parser.add_option('-f', '--fmin', dest='fmin', action='store', default='0.01',
            help='Minimum frequency', metavar='FMIN')
    parser.add_option('-F', '--fmax', dest='fmax', action='store', default='50.0',
            help='Maximum frequency', metavar='FMAX')
    parser.add_option('-k', '--fc', dest='fc', action='store', default='10.0',
            help='Corner frequency', metavar='FC')
    parser.add_option('-m', '--mag', dest='mag', action='store', default='2.0',
            help='Moment magnitude', metavar='Mw')
    parser.add_option('-M', '--moment', dest='Mo', action='store', default='NaN',
            help='Seismic moment', metavar='Mo')
    parser.add_option('-t', '--tstar', dest='t_star', action='store', default='0.0',
            help='T-star (attenuation)', metavar='T')
    parser.add_option('-C', '--combine', dest='combine', action='store_true', default=False,
            help='Generate all the combinations of fc, mag, Mo, tstar')
    parser.add_option('-p', '--plot', dest='plot', action='store_true', default=False,
            help='Plot results')


    (options, args) = parser.parse_args()
    if len(args) < 0:
        parser.print_usage(file=sys.stderr)
        sys.stderr.write("\tUse '-h' for help\n\n")
        sys.exit(1)

    options.fmin = float(options.fmin)
    options.fmax = float(options.fmax)

    options.mag = map(float, options.mag.rstrip(',').split(','))
    options.Mo = map(float, options.Mo.rstrip(',').split(','))
    #A#{
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
        oplist = [(fc, mag, Mo, t_star) for fc in options.fc for mag in options.mag for Mo in options.Mo for t_star in options.t_star]
        oplist = map(list, zip(*oplist))
    else:
        # Add trailing "None" to shorter lists and zip:
        oplist = [options.fc, options.mag, options.Mo, options.t_star]
        oplist = map(None, *oplist)
        # Unzip and convert tuple to lists:
        oplist = map(list, zip(*oplist))
        for l in oplist:
            for n, x in enumerate(l):
                if x == None:
                    l[n] = l[n-1]

    options.fc, options.mag, options.Mo, options.t_star = oplist
    #A#}
    # Add unused options (required by source_spec):
    options.pick_file = None

    return options, args

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
    c.write(open(configfile, 'w'))
    print 'Sample config file written to: ' + configfile

def configure(progname='source_spec'):
    global DEBUG
    if progname == 'source_spec':
        options, args = __parse_args_source_spec()
    elif progname == 'source_model':
        options, args = __parse_args_source_model()
    else:
        sys.stderr.write('Wrong program name: %s\n' % progname)
        sys.exit(1)

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

    #add options and args to config:
    config.options = options
    config.args = args

    return config

oldlogfile=None
def setup_logging(config, basename=None):
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
