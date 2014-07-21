# -*- coding: utf8 -*-
# ssp_residuals.py
#
# (c) 2013-2014 Claudio Satriano <satriano@ipgp.fr>,
#               Agnes Chounet <chounet@ipgp.fr>
'''
Spectral residual routine for source_spec.
'''
import os
import logging
import cPickle as pickle
from obspy.core import Stream
from ssp_spectral_model import spectral_model
from ssp_util import mag_to_moment

def spectral_residuals(config, spec_st, evid, sourcepar_mean):
    '''
    Computes spectral residuals with respect to an average spectral
    model and saves a stream of residuals to disk using pickle.
    '''

    residuals = Stream()
    for station in set(x.stats.station for x in spec_st.traces):
        spec_st_sel = spec_st.select(station=station)
        for spec in spec_st_sel.traces:
            if spec.stats.channel[-1] != 'H':
                continue

            xdata = spec.get_freq()
            synth_mean_mag = spectral_model(xdata, **sourcepar_mean)

            res = spec.copy()
            res.data_mag = spec.data_mag - synth_mean_mag
            res.data = mag_to_moment(res.data_mag)
            residuals.append(res)

    # Save residuals as pickle file
    res_file = os.path.join(config.options.outdir, evid + '-residuals.pickle')
    logging.info('Spectral residuals saved to: %s' % res_file)
    with open(res_file, 'wb') as fp:
        pickle.dump(residuals, fp)
