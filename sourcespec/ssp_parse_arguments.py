# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Argument parser for sourcespec.

:copyright:
    2021-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from itertools import zip_longest
from sourcespec._version import get_versions


def _parse_values(value_str):
    # Lazy-import for speed
    import numpy as np
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
            sys.stderr.write('ERROR: Invalid value: {}\n'.format(value_str))
            sys.exit(1)
    return output


def _get_description(progname):
    if progname == 'source_spec':
        nargs = '+'
        description = 'Estimation of seismic source parameters from '
        description += 'inversion of P- or S-wave displacement spectra.'
        epilog = 'Check the online documentation at ' +\
            'https://sourcespec.rtfd.io for more details.'
    elif progname == 'source_model':
        nargs = '*'
        description = 'Direct modeling of P- or S-wave displacement spectra.'
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
        sys.stderr.write('Wrong program name: {}\n'.format(progname))
        sys.exit(1)
    return description, epilog, nargs


def _init_parser(description, epilog, nargs):
    parser = ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=RawTextHelpFormatter
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '-S', '--sampleconf', dest='sampleconf',
        action='store_true', default=False,
        help='write sample configuration to file and exit'
    )
    group.add_argument(
        '-U', '--updateconf', dest='updateconf',
        action='store', default=None,
        help='update an existing config file from a previous version',
        metavar='FILE'
    )
    parser.add_argument(
        '-c', '--configfile', dest='config_file',
        action='store', default='source_spec.conf',
        help='load configuration from FILE (default: source_spec.conf)',
        metavar='FILE'
    )
    group.add_argument(
        '-t', '--trace_path', nargs=nargs,
        help='path to trace file(s) or trace dir'
    )
    parser.add_argument(
        '-q', '--qmlfile', dest='qml_file',
        action='store', default=None,
        help='get picks and hypocenter information from QuakeML FILE',
        metavar='FILE'
    )
    parser.add_argument(
        '-H', '--hypocenter', dest='hypo_file',
        action='store', default=None,
        help='get hypocenter information from FILE.\n'
             'Supported formats: HYPO71, HYPOINVERSE-2000.\n'
             'If this file contains picks, they will be parsed as well.',
        metavar='FILE'
    )
    parser.add_argument(
        '-p', '--pickfile', dest='pick_file',
        action='store', default=None,
        help='get picks from FILE. Supported formats: HYPO71',
        metavar='FILE'
    )
    parser.add_argument(
        '-n', '--evname', dest='evname',
        action='store', default=None,
        help='event name (used for plots and output files) ',
        metavar='NAME'
    )
    parser.add_argument(
        '-e', '--evid', dest='evid',
        action='store', default=None,
        help='get evid from catalog', metavar='EVID'
    )
    parser.add_argument(
        '-s', '--station', dest='station',
        action='store', default=None,
        help='only use this station', metavar='STATION'
    )
    parser.add_argument(
        '-N', '--no-response', dest='no_response',
        action='store_true', default=False,
        help='do not remove instrument response'
    )
    return parser


def _update_parser(parser, progname):
    if progname == 'source_spec':
        parser.add_argument(
            '-o', '--outdir', dest='outdir',
            action='store', default='sspec_out',
            help='save output to OUTDIR (default: sspec_out)',
            metavar='OUTDIR'
        )
        parser.add_argument(
            '-r', '--run_id', dest='run_id',
            action='store', default='',
            help='string identifying current run (default: "")',
            metavar='RUN_ID'
        )
    elif progname == 'source_model':
        parser.add_argument(
            '-f', '--fmin', dest='fmin', action='store',
            type=float, default='0.01',
            help='minimum frequency (Hz, default 0.01)',
            metavar='FMIN'
        )
        parser.add_argument(
            '-F', '--fmax', dest='fmax', action='store',
            type=float, default='50.0',
            help='maximum frequency (Hz, default 50.0)',
            metavar='FMAX'
        )
        parser.add_argument(
            '-k', '--fc', dest='fc', action='store',
            default='10.0',
            help='(list of) corner frequency (Hz, default 10.0)',
            metavar='Fc'
        )
        parser.add_argument(
            '-m', '--mag', dest='mag', action='store',
            default='2.0',
            help='(list of) moment magnitude (default 2.0)',
            metavar='Mw'
        )
        parser.add_argument(
            '-M', '--moment', dest='Mo', action='store',
            default='NaN',
            help='(list of) seismic moment (N.m, default undefined)',
            metavar='Mo'
        )
        parser.add_argument(
            '-*', '--tstar', dest='t_star', action='store',
            default='0.0',
            help='(list of) t-star (attenuation, default 0.0)',
            metavar='T-STAR'
        )
        parser.add_argument(
            '-a', '--alpha', dest='alpha', action='store',
            default='1.0',
            help='(list of) alpha (exponent for frequency dependence\n'
                 'of attenuation, default 1.0)',
            metavar='1.0'
        )
        parser.add_argument(
            '-C', '--combine', dest='combine',
            action='store_true', default=False,
            help='generate all the combinations of fc, mag, Mo, tstar, alpha'
        )
        parser.add_argument(
            '-P', '--plot', dest='plot', action='store_true',
            default=False, help='plot results'
        )
    parser.add_argument(
        '-v', '--version', action='version',
        version=get_versions()['version']
    )


def parse_args(progname):
    """Parse command line arguments."""
    description, epilog, nargs = _get_description(progname)
    parser = _init_parser(description, epilog, nargs)
    _update_parser(parser, progname)

    options = parser.parse_args()

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
            for opt in oplist:
                for n, x in enumerate(opt):
                    if x is None:
                        opt[n] = opt[n-1]

        options.fc, options.mag, options.Mo, options.t_star, options.alpha = \
            oplist
        # Add unused options (required by source_spec):
        options.correction = False
        options.outdir = '.'

    return options
