# -*- coding: utf-8 -*-
# ssp_output.py
#
# Output functions for source_spec
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
import os
import logging
import numpy as np
from ssp_setup import ssp_exit

def write_output(config, evid, sourcepar):
    if len(sourcepar) == 0:
        logging.info('No source parameter calculated') 
        ssp_exit()

    # Write results to the output dir
    if not os.path.exists(config.options.outdir):
        os.makedirs(config.options.outdir)
    parfilename = '%s/%s.ssp.out' % (config.options.outdir, evid)
    parfile = open(parfilename, 'w')

    # Write source parameters
    parfile.write('*** Source parameters ***\n')
    for statId in sorted(sourcepar.keys()):
        par = sourcepar[statId]
        parfile.write('%s\t' % statId)
        for key in par:
            parfile.write('  %s %6.3f ' % (key, par[key]))
        parfile.write('\n')
    parfile.write('\n')

    ## Compute and write average source parameters
    parfile.write('*** Average source parameters ***\n')
    # Mw
    Mw_array = np.array(list(x['Mw'] for x in sourcepar.values()))
    Mw_mean  = Mw_array.mean()
    Mw_std   = Mw_array.std()
    parfile.write('Mw: %.2f +/- %.2f\n' % (Mw_mean, Mw_std))

    # Mo (N.m)
    Mo_array = np.power(10, Mw_array*1.5 + 9.1)
    Mo_mean  = Mo_array.mean()
    Mo_std   = Mo_array.std()
    parfile.write('Mo: %.3e +/- %.3e N.m\n' % (Mo_mean, Mo_std))

    # fc , hertz
    fc_array = np.array(list(x['fc'] for x in sourcepar.values()))
    fc_mean  = fc_array.mean()
    fc_std   = fc_array.std()
    parfile.write('fc: %.3f +/- %.3f Hz\n' % (fc_mean, fc_std))

    # ra, radius (meters)
    vs_m = config.vs*1000
    ra_array = 0.37 * vs_m / fc_array
    ra_mean  = ra_array.mean()
    ra_std   = ra_array.std()
    parfile.write('Source radius: %.3f +/- %.3f m\n' % (ra_mean, ra_std))

    # bds, Brune stress drop (MPa)
    bsd_array = 0.4375 * Mo_array / np.power(ra_array, 3) * 1e-6
    bsd_mean  = bsd_array.mean()
    bsd_std   = bsd_array.std()
    parfile.write('Brune stress drop: %.3f +/- %.3f MPa\n' % (bsd_mean, bsd_std))

    parfile.close()
    logging.info('Output written to: ' + parfilename)


    if config.options.hypo_file:
        fp = open(config.options.hypo_file, 'r')
        line = fp.readline()
        # Check if first 10 digits of the line contain characters
        if any(c.isalpha() for c in line[0:10]):
            line1 = line
            line = fp.readline()
        line = list(line)
        fp.close()
        mag='%03.2f' % Mw_mean
        for i in range(0,4):
            #line[49+i] = mag[0+i]
            line[45+i] = mag[0+i]
        outline = ''.join(line)
        hypo_file_out = '%s/%s.ssp.h' % (config.options.outdir, evid)
        fp = open(hypo_file_out, 'w')
        try: fp.write(line1)
        except: pass
        fp.write(outline)
        fp.close()
        logging.info('Hypo file written to: ' + hypo_file_out)
