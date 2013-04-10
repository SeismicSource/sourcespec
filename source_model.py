#!/usr/bin/env python
# -*- coding: utf8 -*- 
from optparse import OptionParser
from obspy.core import Stream
from lib.ssp_spectral_model import spectral_model
from lib.ssp_plot_spectra import *
from lib.spectrum import Spectrum

def __parse_args():
    usage = 'usage: %prog [options] trace_file(s) | trace_dir'

    parser = OptionParser(usage=usage);
    parser.add_option('-f', '--fmin', dest='fmin', action='store', default='0.01',
            help='Minimum frequency', metavar='FMIN')
    parser.add_option('-F', '--fmax', dest='fmax', action='store', default='50.0',
            help='Maximum frequency', metavar='FMAX')
    parser.add_option('-c', '--fc', dest='fc', action='store', default='10.0',
            help='Corner frequency', metavar='FC')
    parser.add_option('-m', '--mag', dest='mag', action='store', default='2.0',
            help='Moment magnitude', metavar='Mw')
    parser.add_option('-t', '--tstar', dest='t_star', action='store', default='0.0',
            help='T-star (attenuation)', metavar='T')

    (options, args) = parser.parse_args()
    if len(args) < 0:
        parser.print_usage(file=sys.stderr)
        sys.stderr.write("\tUse '-h' for help\n\n")
        sys.exit(1)

    return options, args

def main():
    (options, args) = __parse_args()
    fdelta = 0.01
    fmin = float(options.fmin)
    fmax = float(options.fmax) + fdelta
    fc = float(options.fc)
    mag = float(options.mag)
    t_star = float(options.t_star)
    
    spec = Spectrum()
    spec.stats.begin = fmin
    spec.stats.delta = fdelta
    spec.stats.npts = int((fmax-fmin)/fdelta)
    spec.stats.instrtype = 'Synth'
    spec.stats.channel = 'Synth'
    
    freq = spec.get_freq()
    spec.data = spectral_model(freq, mag, fc, t_star)
    
    spec_st = Stream()
    spec_st.append(spec)
    
    class C(): pass
    config = C()
    config.PLOT_SHOW = True
    config.PLOT_SAVE = False
    
    plot_spectra(config, spec_st, ncols=1)

if __name__ == '__main__':
    main()
