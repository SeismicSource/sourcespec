# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
plot_sourcepars.py

1D or 2D plot of source parameters from a sqlite parameter file.

:copyright:
    2023-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import argparse
import contextlib
import itertools
import sys
import os
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
import scipy.stats as stats
from scipy.optimize import curve_fit
from obspy import UTCDateTime

valid_plot_types = [
    'fc', 'Er', 'ssd', 'ra', 'Mo', 't_star', 'Qo', 'sigma_a',
    'fc_mw', 'Er_mw', 'ssd_mw', 'ssd_depth', 'mw_time', 'GR'
]
# the parameter base color is the same as in ssp_plot_params_stats.py
valid_parameters = {
    'depth': ('Depth', 'km', 'lin', '#1E90FF'),
    'nsta': ('Number of Stations', None, 'lin', '#F0C300'),
    'Mo': ('Seismic Moment', 'N·m', 'log', '#EE5835'),
    'Mw': ('Moment Magnitude', None, 'lin', '#EE5835'),
    'fc': ('Corner Frequency', 'Hz', 'log', '#6FBA6C'),
    'ssd': ('Static Stress Drop', 'MPa', 'log', '#D4ADD2'),
    'ra': ('Source Radius', 'm', 'lin', '#FAAC64'),
    't_star': ('T star', 's', 'lin', '#9EBAE2'),
    'Qo': ('Qo', None, 'lin', '#C07131'),
    'Er': ('Radiated Energy', 'N·m', 'log', '#00E3E9'),
    'sigma_a': ('Apparent Stress', 'MPa', 'log', '#943B99'),
}


def get_param_label(parameter):
    """
    Get the label for the parameter.
    """
    p_name, p_units = valid_parameters[parameter][:2]
    return f'{p_name} ({p_units})' if p_units is not None else p_name


def get_colormap(parameter):
    """
    Get a colormap for the parameter.
    """
    _, _, scale, base_color = valid_parameters[parameter]
    cmap = plt.cm.colors.LinearSegmentedColormap.from_list(
        'custom', ['black', base_color, 'white'])
    return cmap, scale


class Annot():
    """
    Annotate the plot with the evid, Mw and the value of the parameter.
    """
    def __init__(self, xdata, ydata, labels, yformat, xformat=None):
        self.xdata = xdata
        self.ydata = ydata
        self.labels = labels
        self.yformat = yformat
        self.xformat = xformat

    def __call__(self, event):
        ind = event.ind
        x = np.take(self.xdata, ind)
        y = np.take(self.ydata, ind)
        label = np.take(self.labels, ind)
        for xvalue, yvalue, evid in zip(x, y, label):
            xstring = (
                self.xformat.format(xvalue)
                if self.xformat is not None
                else f'Mw {xvalue:.1f}'
            )
            ystring = self.yformat.format(yvalue)
            print(f'{evid} {xstring} {ystring}')


def mag_to_moment(mag, b=0.5):
    """
    Convert magnitude to moment.

    The parameter b is used to change the slope of the fc-Mw stress drop curve.

    The standard value of b is 0.5, which corresponds to a self-similar model.
    """
    mag = np.asarray(mag)
    return np.power(10, (3 * b * mag + 9.1))


def moment_to_mag(moment, b=0.5):
    """
    Convert moment to magnitude.

    The parameter b is used to change the slope of the fc-Mw stress drop curve.

    The standard value of b is 0.5, which corresponds to a self-similar model.
    """
    moment = np.asarray(moment)
    return (np.log10(moment) - 9.1) / (3 * b)


def stress_drop_curve_fc_mw(delta_sigma, vel, mw, k=0.3724, b=-0.5):
    """
    Constant stress drop curve in fc vs Mw.

    Obtained by combining the equation for stress drop:

        delta_sigma = 7/16 * Mo / a^3

    with the equation for source radius:

        a = k * vel * (delta_sigma / fc)

    where k is a coefficient discussed in Kaneko and Shearer (2014).

    For the Brune source model, k=0.3724.

    Parameters
    ----------
    delta_sigma : float
        Stress drop in MPa.
    vel : float
        P or S-wave velocity in km/s.
    mw : float
        Moment magnitude.
    b : float, optional
        The slope of the stress drop curve. Default is -0.5.
    k : float, optional
        Coefficient for the source radius. Default is 0.3724
        (Brune source model).

    Returns
    -------
    fc : float
        Corner frequency in Hz.
    """
    vel *= 1e3  # P or S-wave velocity in m/s
    delta_sigma *= 1e6  # stress drop in Pa
    # compute moment from magnitude, allowing for non-self-similarity through
    # the b parameter
    moment = mag_to_moment(mw, -b)
    # return fc in Hz
    return k*vel*(16/7 * delta_sigma / moment)**(1/3)


def stress_drop_curve_Er_mw(delta_sigma, mu, mw):
    """
    Constant stress drop curve in Er vs Mw.

    Madariaga (2009), doi:10.1007/978-1-4419-7695-6_22, eq. 33., page 374.
    """
    # Eq. (33) of Madariaga (2009):
    #   0.2331 * delta_sigma = mu * Er / Mo
    Mo = mag_to_moment(mw)
    delta_sigma *= 1e6
    # return Er in N·m
    return 0.2331 * delta_sigma * Mo / mu


def apparent_stress_curve_Er_mw(sigma_a, mu, mw):
    """
    Constant apparent stress curve in Er vs Mw.

    Madariaga (2009), doi:10.1007/978-1-4419-7695-6_22, eq. 33., page 374.
    """
    # Eq. (33) of Madariaga (2009):
    #   sigma_a = mu * Er / Mo
    Mo = mag_to_moment(mw)
    sigma_a *= 1e6
    # return Er in N·m
    return sigma_a * Mo / mu


def fc_mw_function(mw, a, b):
    """Function to fit fc vs Mw."""
    return a / 3. + b * mw


def calc_r2(x, y, yerr, a, b):
    """Coefficient of determination."""
    y_mean = np.mean(y)
    y_calc = fc_mw_function(x, a, b)
    SS_tot = np.sum(((y - y_mean) / yerr)**2.)
    SS_res = np.sum(((y - y_calc) / yerr)**2.)
    return 1 - SS_res / SS_tot


def compute_mc(magnitudes, magn_bin, vector_compl, pval):
    """
    Compute the magnitude of completeness of a seismic catalog using the
    method in Taroni 2023 - TSR, https://doi.org/10.1785/0320230017.

    This code is translated from the original MATLAB code provided by
    Matteo Taroni.

    Parameters
    ----------
    magnitudes : np.ndarray
        Vector of magnitudes.
    magn_bin : float
        Binning of the magnitudes (e.g., 0.01 or 0.1).
    vector_compl : np.ndarray
        Vector of completeness values to analyze.
    pval : float
        P-value threshold for the complete part of the catalog (e.g., 0.1).

    Returns
    -------
    mc : float
        The magnitude of completeness of the catalog.
    """
    # Preallocate P-value matrix
    pvalue_unif = np.zeros((100, len(vector_compl)))
    # CDF of uniform distribution
    cdf_unif = stats.uniform()
    # Loop over different completeness values
    for i, j in itertools.product(range(len(vector_compl)), range(100)):
        # Add uniform error to magnitudes
        magn = magnitudes + (np.random.rand(len(magnitudes)) - 0.5) * magn_bin
        # Select magnitudes above completeness threshold
        magn_ok = magn[magn >= vector_compl[i]] - vector_compl[i]
        # Transformation from exponential to uniform random variables
        if len(magn_ok) < 2:
            continue
        transf = magn_ok[:-1:2] / (magn_ok[:-1:2] + magn_ok[1::2])
        # Kolmogorov-Smirnov test for uniformity
        _, pvalue_unif[j, i] = stats.kstest(transf, cdf_unif.cdf)
    p_mean = np.mean(pvalue_unif, axis=0)
    # Find the magnitude of completeness (first mean p-value >= Pval)
    mc_candidates = vector_compl[p_mean >= pval]
    return np.min(mc_candidates) if mc_candidates.size > 0 else None


def fit_gutenberg_richter(bin_centers, cum_nevs, mc):
    """
    Fit the Gutenberg-Richter law to the data points.

    Parameters
    ----------
    bin_centers : np.ndarray
        Bin centers.
    cum_nevs : np.ndarray
        Cumulative number of events.
    mc : float
        Magnitude of completeness.

    Returns
    -------
    a, b : float
        G-R law parameters.
    """
    def gr_law(mw, a, b):
        return a - b * mw
    # discard bins with no events
    log_cum_nevs = np.log10(cum_nevs[cum_nevs > 0])
    bin_centers_fit = bin_centers[cum_nevs > 0]
    # discard bins below the magnitude of completeness
    log_cum_nevs_fit = log_cum_nevs[bin_centers_fit >= mc]
    bin_centers_fit = bin_centers_fit[bin_centers_fit >= mc]
    # with full_output = False, always returns a 2-tuple
    # pylint: disable-next=unbalanced-tuple-unpacking
    popt, _ = curve_fit(
        gr_law, bin_centers_fit, log_cum_nevs_fit,
        full_output=False
    )
    a, b = popt
    # discard bins which are too far away from the fit and re-fit
    while True:
        residuals = log_cum_nevs_fit - gr_law(bin_centers_fit, a, b)
        std_res = np.std(residuals)
        cond = np.abs(residuals) < 2 * std_res
        if np.all(cond):
            break
        log_cum_nevs_fit = log_cum_nevs_fit[cond]
        bin_centers_fit = bin_centers_fit[cond]
        # pylint: disable-next=unbalanced-tuple-unpacking
        popt, _ = curve_fit(
            gr_law, bin_centers_fit, log_cum_nevs_fit, p0=(a, b),
            full_output=False
        )
        a, b = popt
    return a, b


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.description =\
        '1D or 2D plot of source parameters from a sqlite parameter file'
    parser.add_argument('sqlite_file')
    valid_plot_types_str = str(valid_plot_types)[1:-1].replace("'", "\"")
    parser.add_argument(
        '-p', '--plot_type', default='fc_mw',
        help=f'Plot type. One of: {valid_plot_types_str}. '
             '1D plots are histograms. '
             '2D plots are scatter plots or 2D histograms, '
             'depending on the -H option. '
             'Default is "fc_mw"'
    )
    valid_parameters_str = str(list(
        valid_parameters.keys()))[1:-1].replace("'", "\"")
    parser.add_argument(
        '-c', '--colorby', default=None,
        help='Color the data points by this parameter. '
             f'One of: {valid_parameters_str}. '
             'Default is None')
    parser.add_argument(
        '-C', '--colormap', default=None,
        help='Force using a specific colormap instead of the default one for '
             'the colorby parameter. Any Matplotlib colormap is accepted. ')
    parser.add_argument(
        '-r', '--runid', default=None,
        help='Only select a specific runid. Default is all. Use "latest" to '
             'select the latest runid for each event')
    parser.add_argument(
        '-s', '--statistics', default='mean',
        help='Statistics to use: "mean", "wmean" (weighted mean) '
             'or "pctl" (percentiles). Default is "mean"')
    parser.add_argument(
        '-i', '--nbins', type=int, default=None,
        help='Number of bins in the histogram (default: autoset)')
    parser.add_argument(
        '-n', '--stamin', type=int, default=None,
        help='Minimum number of stations')
    parser.add_argument(
        '-N', '--stamax', type=int, default=None,
        help='Maximum number of stations')
    parser.add_argument(
        '-m', '--magmin', type=float, default=None,
        help='Minimum magnitude')
    parser.add_argument(
        '-M', '--magmax', type=float, default=None,
        help='Maximum magnitude')
    parser.add_argument(
        '-d', '--ssdmin', type=float, default=None,
        help='Minimum static stress drop')
    parser.add_argument(
        '-D', '--ssdmax', type=float, default=None,
        help='Maximum static stress drop')
    parser.add_argument(
        '-a', '--sigmaamin', type=float, default=None,
        help='Minimum apparent stress drop')
    parser.add_argument(
        '-A', '--sigmaamax', type=float, default=None,
        help='Maximum apparent stress drop')
    parser.add_argument(
        '-H', '--hist', default=False, action='store_true',
        help='Draw an histogram instead of a scatter plot (only for 2D plots)')
    parser.add_argument(
        '-f', '--fit', default=False, action='store_true',
        help='Fit the scatter plot with a linear function (only for 2D plots)')
    parser.add_argument(
        '-l', '--slope', default=False, action='store_true',
        help='Fit also the slope of the linear function (only for 2D plots)')
    parser.add_argument(
        '-w', '--wave_type', default='S',
        help='Wave type. Only used for plots involving "fc". '
             'One of "P", "S", "SV" or "SH". Default is "S".')
    args = parser.parse_args()
    if args.plot_type not in valid_plot_types:
        msg = f'Plot type must be one of {valid_plot_types_str}'
        raise ValueError(msg)
    if args.colorby is not None and args.colorby not in valid_parameters:
        msg = f'Colorby parameter must be one of {valid_parameters_str}'
        raise ValueError(msg)
    if args.statistics not in ['mean', 'wmean', 'pctl']:
        msg = 'Statistics must be "mean", "wmean" or "pctl"'
        raise ValueError(msg)
    if args.wave_type not in ['P', 'S', 'SV', 'SH']:
        msg = 'Wave type must be "P", "S", "SV" or "SH"'
        raise ValueError(msg)
    return args


def query_event_params_into_numpy(cursor, param, param_type, query_condition):
    """
    Query a parameter from the Events table and return it as a numpy array.
    """
    query = (
        f'SELECT e.{param} FROM Events e {query_condition} '
        'ORDER BY e.evid, e.runid'
    )
    cursor.execute(query)
    result = np.array(cursor.fetchall())
    if len(result) == 0:
        raise ValueError('No events found')
    if param_type == UTCDateTime:
        return np.array([UTCDateTime(t) for t in result[:, 0]])
    return result[:, 0].astype(param_type)


class Params():
    """
    Class to handle the parameters from a sqlite file.
    """
    def __init__(self, args):
        """
        Initialize the class from a sqlite file.
        """
        self.args = args
        self.sqlite_file = args.sqlite_file
        self.runid = runid = args.runid
        self.stat = stat = args.statistics
        self._open_db(self.sqlite_file)
        if runid == 'latest':
            # select the latest runid for each event
            query_condition = """
JOIN (
    SELECT evid, MAX(runid) AS max_runid
    FROM events
    GROUP BY evid
) AS max_runids
ON e.evid = max_runids.evid AND e.runid = max_runids.max_runid
"""
        elif runid is not None:
            query_condition = f'WHERE runid="{runid}"'
        else:
            query_condition = ''
        self.evids = query_event_params_into_numpy(
            self.cur, 'evid', str, query_condition)
        self.time = query_event_params_into_numpy(
            self.cur, 'orig_time', UTCDateTime, query_condition)
        self.depth = query_event_params_into_numpy(
            self.cur, 'depth', np.float64, query_condition)
        self.vp = query_event_params_into_numpy(
            self.cur, 'vp', np.float64, query_condition)
        self.vs = query_event_params_into_numpy(
            self.cur, 'vs', np.float64, query_condition)
        self.rho = query_event_params_into_numpy(
            self.cur, 'rho', np.float64, query_condition)
        self.kp = query_event_params_into_numpy(
            self.cur, 'kp', np.float64, query_condition)
        self.ks = query_event_params_into_numpy(
            self.cur, 'ks', np.float64, query_condition)
        self.wave_type = query_event_params_into_numpy(
            self.cur, 'wave_type', str, query_condition)
        self.nsta = query_event_params_into_numpy(
            self.cur, 'nobs', np.int32, query_condition)
        self.Mo = query_event_params_into_numpy(
            self.cur, f'Mo_{stat}', np.float64, query_condition)
        self.Mw = query_event_params_into_numpy(
            self.cur, f'Mw_{stat}', np.float64, query_condition)
        self.Mw_err_minus = query_event_params_into_numpy(
            self.cur, f'Mw_{stat}_err_minus', np.float64, query_condition)
        self.Mw_err_plus = query_event_params_into_numpy(
            self.cur, f'Mw_{stat}_err_plus', np.float64, query_condition)
        self.fc = query_event_params_into_numpy(
            self.cur, f'fc_{stat}', np.float64, query_condition)
        self.fc_err_minus = query_event_params_into_numpy(
            self.cur, f'fc_{stat}_err_minus', np.float64, query_condition)
        self.fc_err_plus = query_event_params_into_numpy(
            self.cur, f'fc_{stat}_err_plus', np.float64, query_condition)
        self.Er = query_event_params_into_numpy(
            self.cur, f'Er_{stat}', np.float64, query_condition)
        self.Er_err_minus = query_event_params_into_numpy(
            self.cur, f'Er_{stat}_err_minus', np.float64, query_condition)
        self.Er_err_plus = query_event_params_into_numpy(
            self.cur, f'Er_{stat}_err_plus', np.float64, query_condition)
        self.ssd = query_event_params_into_numpy(
            self.cur, f'ssd_{stat}', np.float64, query_condition)
        self.ssd_err_minus = query_event_params_into_numpy(
            self.cur, f'ssd_{stat}_err_minus', np.float64, query_condition)
        self.ssd_err_plus = query_event_params_into_numpy(
            self.cur, f'ssd_{stat}_err_plus', np.float64, query_condition)
        self.ra = query_event_params_into_numpy(
            self.cur, f'ra_{stat}', np.float64, query_condition)
        self.ra_err_minus = query_event_params_into_numpy(
            self.cur, f'ra_{stat}_err_minus', np.float64, query_condition)
        self.ra_err_plus = query_event_params_into_numpy(
            self.cur, f'ra_{stat}_err_plus', np.float64, query_condition)
        self.sigma_a = query_event_params_into_numpy(
            self.cur, f'sigma_a_{stat}', np.float64, query_condition)
        self.sigma_a_err_minus = query_event_params_into_numpy(
            self.cur, f'sigma_a_{stat}_err_minus', np.float64, query_condition)
        self.sigma_a_err_plus = query_event_params_into_numpy(
            self.cur, f'sigma_a_{stat}_err_plus', np.float64, query_condition)
        self.t_star = query_event_params_into_numpy(
            self.cur, f't_star_{stat}', np.float64, query_condition)
        self.t_star_err_minus = query_event_params_into_numpy(
            self.cur, f't_star_{stat}_err_minus', np.float64, query_condition)
        self.t_star_err_plus = query_event_params_into_numpy(
            self.cur, f't_star_{stat}_err_plus', np.float64, query_condition)
        self.Qo = query_event_params_into_numpy(
            self.cur, f'Qo_{stat}', np.float64, query_condition)
        self.Qo_err_minus = query_event_params_into_numpy(
            self.cur, f'Qo_{stat}_err_minus', np.float64, query_condition)
        self.Qo_err_plus = query_event_params_into_numpy(
            self.cur, f'Qo_{stat}_err_plus', np.float64, query_condition)
        self.cur.close()
        # other attributes
        self.nbins_x = None
        self.nbins_y = None

    def _open_db(self, sqlite_file):
        """
        Open the sqlite file and check that it contains the required tables.
        """
        if not os.path.isfile(sqlite_file):
            raise FileNotFoundError(f'File "{sqlite_file}" not found')
        if os.stat(sqlite_file).st_size == 0:
            raise ValueError(f'File "{sqlite_file}" is empty')
        conn = sqlite3.connect(sqlite_file)
        cur = conn.cursor()
        try:
            cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
        except sqlite3.DatabaseError as e:
            raise ValueError(
                f'File "{sqlite_file}" is not a valid sqlite file'
            ) from e
        tables = [t[0] for t in cur.fetchall()]
        for table in 'Events', 'Stations':
            if table not in tables:
                raise ValueError(
                    f'Table "{table}" not found in file "{sqlite_file}"')
        self.cur = cur

    def skip_events(self, idx):
        """Skip events with index idx."""
        self.evids = np.delete(self.evids, idx)
        self.time = np.delete(self.time, idx)
        self.depth = np.delete(self.depth, idx)
        self.vp = np.delete(self.vp, idx)
        self.vs = np.delete(self.vs, idx)
        self.rho = np.delete(self.rho, idx)
        self.kp = np.delete(self.kp, idx)
        self.ks = np.delete(self.ks, idx)
        self.wave_type = np.delete(self.wave_type, idx)
        self.nsta = np.delete(self.nsta, idx)
        self.Mo = np.delete(self.Mo, idx)
        self.Mw = np.delete(self.Mw, idx)
        self.Mw_err_minus = np.delete(self.Mw_err_minus, idx)
        self.Mw_err_plus = np.delete(self.Mw_err_plus, idx)
        self.fc = np.delete(self.fc, idx)
        self.fc_err_minus = np.delete(self.fc_err_minus, idx)
        self.fc_err_plus = np.delete(self.fc_err_plus, idx)
        self.Er = np.delete(self.Er, idx)
        self.Er_err_minus = np.delete(self.Er_err_minus, idx)
        self.Er_err_plus = np.delete(self.Er_err_plus, idx)
        self.ssd = np.delete(self.ssd, idx)
        self.ssd_err_minus = np.delete(self.ssd_err_minus, idx)
        self.ssd_err_plus = np.delete(self.ssd_err_plus, idx)
        self.ra = np.delete(self.ra, idx)
        self.ra_err_minus = np.delete(self.ra_err_minus, idx)
        self.ra_err_plus = np.delete(self.ra_err_plus, idx)
        self.sigma_a = np.delete(self.sigma_a, idx)
        self.sigma_a_err_minus = np.delete(self.sigma_a_err_minus, idx)
        self.sigma_a_err_plus = np.delete(self.sigma_a_err_plus, idx)
        self.t_star = np.delete(self.t_star, idx)
        self.t_star_err_minus = np.delete(self.t_star_err_minus, idx)
        self.t_star_err_plus = np.delete(self.t_star_err_plus, idx)
        self.Qo = np.delete(self.Qo, idx)
        self.Qo_err_minus = np.delete(self.Qo_err_minus, idx)
        self.Qo_err_plus = np.delete(self.Qo_err_plus, idx)

    def filter(self, stamin=None, stamax=None, magmin=None, magmax=None,
               ssdmin=None, ssdmax=None, sigmaamin=None, sigmaamax=None):
        """Filter the parameters based on one or more conditions."""
        cond = np.ones(len(self.nsta)).astype(bool)
        if stamin is not None:
            cond = np.logical_and(cond, self.nsta >= stamin)
        if stamax is not None:
            cond = np.logical_and(cond, self.nsta <= stamax)
        if magmin is not None:
            cond = np.logical_and(cond, self.Mw >= magmin)
        if magmax is not None:
            cond = np.logical_and(cond, self.Mw <= magmax)
        if ssdmin is not None:
            cond = np.logical_and(cond, self.ssd >= ssdmin)
        if ssdmax is not None:
            cond = np.logical_and(cond, self.ssd <= ssdmax)
        if sigmaamin is not None:
            cond = np.logical_and(cond, self.sigma_a >= sigmaamin)
        if sigmaamax is not None:
            cond = np.logical_and(cond, self.sigma_a <= sigmaamax)
        self.skip_events(np.where(~cond)[0])

    def _make_mw_axis(self):
        """Make the magnitude axis."""
        if self.args.magmin is not None:
            mag_min = self.args.magmin
        else:
            mag_min = np.floor(np.min(self.Mw - self.Mw_err_minus))
        if self.args.magmax is not None:
            mag_max = self.args.magmax
        else:
            mag_max = np.ceil(np.max(self.Mw + self.Mw_err_plus))
        xlim_mag = (mag_min, mag_max)
        fig = plt.figure(figsize=(10, 6))
        ax_Mo = fig.add_subplot(111)
        ax_Mo.set_xscale('log')
        ax_Mo.set_xlim(mag_to_moment(xlim_mag))
        ax = ax_Mo.twiny()
        ax.set_xlim(xlim_mag)
        ax.set_xlabel('Mw')
        ax_Mo.set_xlabel(get_param_label('Mo'))
        return fig, ax, ax_Mo

    def _add_grid(self, ax):
        """Add grid to the plot."""
        ax.grid(
            True, which='major', linestyle='solid',
            color='#999999', zorder=0)
        ax.grid(
            True, which='minor', linestyle='solid',
            color='#DDDDDD', zorder=0)

    def _set_plot_title(self, ax, nevs=None, extra_text=None, align='center'):
        """Set the plot title."""
        if nevs is None:
            nevs = len(self.evids)
        stat_descr = {
            'mean': 'mean',
            'wmean': 'weighted mean',
            'pctl': 'percentiles',
        }
        title = f'{nevs} events'
        if None not in (self.nbins_x, self.nbins_y):
            title += f', {self.nbins_x}x{self.nbins_y} bins'
        title += f' - {stat_descr[self.stat]}'
        if self.runid is not None:
            title += f' - runid: {self.runid}'
        if extra_text is not None:
            title += f'\n{extra_text}'
        align_options = {
            'center': (0.5, 0.95, 'center', 'top'),
            'left': (0.05, 0.95, 'left', 'top'),
            'right': (0.95, 0.95, 'right', 'top')
        }
        try:
            xpos, ypos, ha, va = align_options[align]
        except KeyError as e:
            raise ValueError('Invalid align') from e
        bbox = {
            'facecolor': 'white',
            'edgecolor': 'black',
            'boxstyle': 'round',
            'alpha': 0.7
        }
        ax.set_title(
            title, x=xpos, y=ypos,
            horizontalalignment=ha, verticalalignment=va,
            bbox=bbox
        )

    def _stress_drop_curves_fc_mw(self, vel, k_parameter, ax):
        """Plot stress-drop curves for different delta_sigma."""
        mag_min, mag_max = ax.get_xlim()
        mw_step = 0.1
        mw_test = np.arange(mag_min, mag_max - 2 * mw_step, mw_step)
        fc_min = np.inf
        fc_max = 0.
        for delta_sigma in (0.1, 1., 10., 100.):
            fc_test = stress_drop_curve_fc_mw(
                delta_sigma, vel, mw_test, k_parameter)
            if fc_test.min() < fc_min:
                fc_min = fc_test.min()
            if fc_test.max() > fc_max:
                fc_max = fc_test.max()
            ax.plot(mw_test, fc_test, color='#555555')
            label = str(delta_sigma)
            # Get rid of ".0" in label
            if label.endswith('.0'):
                label = label[:-2]
            label += ' MPa'
            outline = patheffects.withStroke(linewidth=2, foreground='white')
            ax.text(mw_test[-1], fc_test[-1], label, path_effects=[outline])
        ax.set_ylim((fc_min * 0.9, fc_max * 1.1))

    def _stress_drop_curves_Er_mw(self, mu, ax):
        """Plot stress-drop curves for different delta_sigma."""
        mag_min, mag_max = ax.get_xlim()
        mw_step = 0.1
        mw_test = np.arange(mag_min, mag_max - 2 * mw_step, mw_step)
        Er_min = np.inf
        Er_max = 0.
        for delta_sigma in (0.1, 1., 10., 100.):
            Er_test = stress_drop_curve_Er_mw(delta_sigma, mu, mw_test)
            if Er_test.min() < Er_min:
                Er_min = Er_test.min()
            if Er_test.max() > Er_max:
                Er_max = Er_test.max()
            ax.plot(mw_test, Er_test, color='#555555')
            label = str(delta_sigma)
            # Get rid of ".0" in label
            if label.endswith('.0'):
                label = label[:-2]
            label += ' MPa'
            outline = patheffects.withStroke(linewidth=2, foreground='white')
            ax.text(mw_test[-1], Er_test[-1], label, path_effects=[outline])
        ax.set_ylim((Er_min * 0.5, Er_max * 2))

    def _apparent_stress_curves_Er_mw(self, mu, ax):
        """Plot apparent stress curves for different delta_sigma."""
        mag_min, mag_max = ax.get_xlim()
        mw_step = 0.1
        mw_test = np.arange(mag_min, mag_max - 2 * mw_step, mw_step)
        Er_min = np.inf
        Er_max = 0.
        for sigma_a in (0.1, 1., 10., 100.):
            Er_test = apparent_stress_curve_Er_mw(sigma_a, mu, mw_test)
            if Er_test.min() < Er_min:
                Er_min = Er_test.min()
            if Er_test.max() > Er_max:
                Er_max = Er_test.max()
            ax.plot(mw_test, Er_test, color='#555555')
            label = str(sigma_a)
            # Get rid of ".0" in label
            if label.endswith('.0'):
                label = label[:-2]
            label += ' MPa'
            outline = patheffects.withStroke(linewidth=2, foreground='white')
            ax.text(mw_test[-1], Er_test[-1], label, path_effects=[outline])
        ax.set_ylim((Er_min * 0.5, Er_max * 2))

    def _2d_hist(self, fig, ax, x, y, x_bins, y_bins, cbaxes_location):
        """General method to plot 2d histograms."""
        counts, _, _ = np.histogram2d(x, y, bins=(x_bins, y_bins))
        cm = ax.pcolormesh(
            x_bins[:-1], y_bins[:-1], counts.T,
            cmap='magma_r', shading='auto')
        if cbaxes_location == 'left':
            cbaxes = fig.add_axes([0.15, 0.15, 0.02, 0.2])
        elif cbaxes_location == 'right':
            cbaxes = fig.add_axes([0.80, 0.15, 0.02, 0.2])
        else:
            raise ValueError('Invalid cbaxes_location')
        plt.colorbar(cm, cax=cbaxes, orientation='vertical', label='counts')
        return len(y)

    def _2d_hist_fc_mw(self, fig, ax, nbins=None, wave_type='S'):
        """Plot a 2d histogram of fc vs mw."""
        mw_min, mw_max = ax.get_xlim()
        fc_min, fc_max = ax.get_ylim()
        log_fc_min = np.log10(fc_min)
        log_fc_max = np.log10(fc_max)
        if nbins is None:
            mw_nbins = int((mw_max - mw_min) * 10)
            fc_nbins = int((log_fc_max - log_fc_min) * 10)
        else:
            mw_nbins = nbins
            fc_nbins = nbins
        self.nbins_x = mw_nbins
        self.nbins_y = fc_nbins
        mw_bins = np.linspace(mw_min, mw_max + 0.1, mw_nbins)
        fc_bins = 10**np.linspace(log_fc_min, log_fc_max + 0.1, fc_nbins)
        cond = self.wave_type == wave_type
        mw = self.Mw[cond]
        if len(mw) == 0:
            raise ValueError(f'No events found for wave type "{wave_type}"')
        fc = self.fc[cond]
        return self._2d_hist(
            fig, ax, mw, fc, mw_bins, fc_bins, cbaxes_location='left')

    def _2d_hist_Er_mw(self, fig, ax, nbins=None):
        """Plot a 2d histogram of Er vs mw."""
        mw_min, mw_max = ax.get_xlim()
        Er_min, Er_max = ax.get_ylim()
        log_Er_min = np.log10(Er_min)
        log_Er_max = np.log10(Er_max)
        if nbins is None:
            mw_nbins = int((mw_max - mw_min) * 10)
            Er_nbins = int((log_Er_max - log_Er_min) * 5)
        else:
            mw_nbins = nbins
            Er_nbins = nbins
        mw_bins = np.linspace(mw_min, mw_max + 0.1, mw_nbins)
        Er_bins = 10**np.linspace(log_Er_min, log_Er_max + 0.1, Er_nbins)
        self.nbins_x = mw_nbins
        self.nbins_y = Er_nbins
        cond = ~np.isnan(self.Er)
        mw = self.Mw[cond]
        Er = self.Er[cond]
        return self._2d_hist(
            fig, ax, mw, Er, mw_bins, Er_bins, cbaxes_location='right')

    def _2d_hist_ssd_mw(self, fig, ax, nbins=None):
        """Plot a 2d histogram of ssd vs mw."""
        mw_min, mw_max = ax.get_xlim()
        ssd_min, ssd_max = ax.get_ylim()
        log_ssd_min = np.log10(ssd_min)
        log_ssd_max = np.log10(ssd_max)
        if nbins is None:
            mw_nbins = int((mw_max - mw_min) * 10)
            ssd_nbins = int((log_ssd_max - log_ssd_min) * 5)
        else:
            mw_nbins = nbins
            ssd_nbins = nbins
        mw_bins = np.linspace(mw_min, mw_max + 0.1, mw_nbins)
        ssd_bins = 10**np.linspace(log_ssd_min, log_ssd_max + 0.1, ssd_nbins)
        self.nbins_x = mw_nbins
        self.nbins_y = ssd_nbins
        cond = ~np.isnan(self.ssd)
        mw = self.Mw[cond]
        ssd = self.ssd[cond]
        return self._2d_hist(
            fig, ax, mw, ssd, mw_bins, ssd_bins, cbaxes_location='left')

    def _2d_hist_ssd_depth(self, fig, ax, nbins=None):
        """Plot a 2d histogram of ssd vs depth."""
        depth_min, depth_max = ax.get_xlim()
        ssd_min, ssd_max = ax.get_ylim()
        log_ssd_min = np.log10(ssd_min)
        log_ssd_max = np.log10(ssd_max)
        if nbins is None:
            depth_nbins = int((depth_max - depth_min) * 5)
            ssd_nbins = int((log_ssd_max - log_ssd_min) * 5)
        else:
            depth_nbins = nbins
            ssd_nbins = nbins
        depth_bins = np.linspace(depth_min, depth_max + 0.1, depth_nbins)
        ssd_bins = 10**np.linspace(log_ssd_min, log_ssd_max + 0.1, ssd_nbins)
        self.nbins_x = depth_nbins
        self.nbins_y = ssd_nbins
        cond = ~np.isnan(self.ssd)
        depth = self.depth[cond]
        ssd = self.ssd[cond]
        return self._2d_hist(
            fig, ax, depth, ssd, depth_bins, ssd_bins, cbaxes_location='left')

    def _scatter_plot(self, fig, ax, x, y, x_err, y_err, evids, yformat, cond,
                      colorby=None, colormap=None, cbaxes_location='left',
                      xformat=None):
        """
        General method to plot scatter plots with error bars and optional
        color coding.
        """
        if colorby is not None:
            color = getattr(self, colorby)[cond]
            cmap, units = get_colormap(colorby)
            if colormap is not None:
                cmap = plt.get_cmap(colormap)
            label = get_param_label(colorby)
            norm = plt.cm.colors.LogNorm() if units == 'log' else None
            mfc = 'none'
            mec = 'none'
            ecolor = 'black'
            errorbar_alpha = 0.5
            scatter_alpha = 1
        else:
            color = 'none'
            cmap = None
            norm = None
            mfc = ecolor = '#FCBA25'
            mec = 'black'
            errorbar_alpha = 1
            scatter_alpha = 0
        ax.errorbar(
            x, y,
            xerr=x_err,
            yerr=y_err,
            fmt='o', mec=mec, mfc=mfc, ecolor=ecolor,
            alpha=errorbar_alpha)
        sc = ax.scatter(
            x, y, alpha=scatter_alpha, c=color, cmap=cmap, norm=norm,
            edgecolors='black',
            picker=True, zorder=20)
        if colorby is not None:
            if cbaxes_location == 'left':
                cbaxes = fig.add_axes([0.15, 0.15, 0.02, 0.2])
            elif cbaxes_location == 'right':
                cbaxes = fig.add_axes([0.80, 0.15, 0.02, 0.2])
            else:
                raise ValueError('Invalid cbaxes_location')
            if len(label) > 10:
                label = label.replace(' (', '\n(')
            cbaxes.set_zorder(1000)
            plt.colorbar(sc, cax=cbaxes, label=label)
        annot = Annot(x, y, evids, yformat, xformat)
        fig.canvas.mpl_connect('pick_event', annot)
        return len(y)

    def _scatter_fc_mw(self, fig, ax, wave_type='S',
                       colorby=None, colormap=None):
        """Plot the scatter plot of fc vs mw."""
        cond = self.wave_type == wave_type
        evids = self.evids[cond]
        if len(evids) == 0:
            raise ValueError(f'No events found for wave type "{wave_type}"')
        mw = self.Mw[cond]
        mw_err_minus = self.Mw_err_minus[cond]
        mw_err_plus = self.Mw_err_plus[cond]
        fc = self.fc[cond]
        fc_err_minus = self.fc_err_minus[cond]
        fc_err_plus = self.fc_err_plus[cond]
        yformat = 'fc {:.2f} Hz'
        return self._scatter_plot(
            fig, ax, mw, fc,
            [mw_err_minus, mw_err_plus], [fc_err_minus, fc_err_plus],
            evids, yformat, cond,
            colorby, colormap, cbaxes_location='left')

    def _scatter_Er_mw(self, fig, ax, colorby=None, colormap=None):
        """Plot the scatter plot of Er vs mw."""
        cond = ~np.isnan(self.Er)
        mw = self.Mw[cond]
        mw_err_minus = self.Mw_err_minus[cond]
        mw_err_plus = self.Mw_err_plus[cond]
        Er = self.Er[cond]
        Er_err_minus = self.Er_err_minus[cond]
        Er_err_plus = self.Er_err_plus[cond]
        evids = self.evids[cond]
        yformat = 'Er {:.1e} N·m'
        return self._scatter_plot(
            fig, ax, mw, Er,
            [mw_err_minus, mw_err_plus],
            [Er_err_minus, Er_err_plus],
            evids, yformat, cond,
            colorby, colormap, cbaxes_location='right')

    def _scatter_ssd_mw(self, fig, ax, colorby=None, colormap=None):
        """Plot the scatter plot of ssd vs mw."""
        cond = ~np.isnan(self.ssd)
        mw = self.Mw[cond]
        mw_err_minus = self.Mw_err_minus[cond]
        mw_err_plus = self.Mw_err_plus[cond]
        ssd = self.ssd[cond]
        ssd_err_minus = self.ssd_err_minus[cond]
        ssd_err_plus = self.ssd_err_plus[cond]
        evids = self.evids[cond]
        yformat = 'ssd {:.2e} MPa'
        return self._scatter_plot(
            fig, ax, mw, ssd,
            [mw_err_minus, mw_err_plus],
            [ssd_err_minus, ssd_err_plus],
            evids, yformat, cond,
            colorby, colormap)

    def _scatter_ssd_depth(self, fig, ax, colorby=None, colormap=None):
        """Plot the scatter plot of ssd vs depth."""
        cond = ~np.isnan(self.ssd)
        depth = self.depth[cond]
        ssd = self.ssd[cond]
        ssd_err_minus = self.ssd_err_minus[cond]
        ssd_err_plus = self.ssd_err_plus[cond]
        evids = self.evids[cond]
        yformat = 'ssd {:.2e} MPa'
        return self._scatter_plot(
            fig, ax, depth, ssd,
            None, [ssd_err_minus, ssd_err_plus],
            evids, yformat, cond,
            colorby, colormap)

    def _fit_fc_mw(self, vel, k_parameter, ax, slope=False):
        """Plot a linear regression of fc vs mw."""
        mag_min, mag_max = ax.get_xlim()
        mw_step = 0.1
        mw_test = np.arange(mag_min, mag_max - 2 * mw_step, mw_step)
        if slope:
            f = fc_mw_function
        else:
            def f(mw, a):
                return fc_mw_function(mw, a, -0.5)
        y = np.log10(self.fc)
        # pylint: disable=unbalanced-tuple-unpacking
        popt, _pcov = curve_fit(f, self.Mw, y)
        try:
            a, b = popt
        except Exception:
            a, = popt
            b = -0.5
        print(f'a: {a:.1f} b {b:.1f}:')
        slope = (3 / 2) / b
        print(f'slope: {slope:.1f}')
        r2 = calc_r2(self.Mw, y, np.ones_like(y), a, b)
        print('r2:', r2)
        # print('r2_err:', r2_err)
        delta_sigma = 1. / ((vel * 1000.)**3.) * 10**(a + 9.1 + 0.935)
        delta_sigma /= 1e6
        print('delta_sigma', delta_sigma)
        fc_test = stress_drop_curve_fc_mw(
            delta_sigma, vel, mw_test, k_parameter, b=b)
        ax.plot(
            mw_test, fc_test, color='red', zorder=10, label=f'r2: {r2:.2f}')
        ax.legend()
        if not slope:
            ax.text(
                mw_test[-1], fc_test[-1], f'{delta_sigma:.1f} MPa',
                color='red', backgroundcolor='white', zorder=50)

    def plot_fc_mw(self, hist=False, fit=False, slope=False, nbins=None,
                   wave_type='S', colorby=None, colormap=None):
        """
        Plot the logarithm of the corner frequency vs the moment magnitude.

        Parameters
        ----------
        hist : bool
            If True, plot a 2D histogram instead of a scatter plot.
        fit : bool
            If True, plot a linear regression of fc vs mw.
        slope : bool
            If True, also fit the slope of the linear regression.
        nbins : int
            Number of bins for the 2D histogram.
        wave_type : str
            One of "P", "S", "SV" or "SH". Default is "S".
        colorby : str
            Color the data points by the given parameter.
        colormap : str
            Use this colormap for the colorby parameter instead of the default.
        """
        fig, ax, ax_Mo = self._make_mw_axis()
        ax_Mo.set_yscale('log')

        if wave_type == 'P':
            vel = np.nanmean(self.vp)
            k = np.nanmean(self.kp)
        elif wave_type in ('S', 'SV', 'SH'):
            vel = np.nanmean(self.vs)
            k = np.nanmean(self.ks)
        else:
            raise ValueError('Wave type must be "P", "S", "SV" or "SH"')
        self._stress_drop_curves_fc_mw(vel, k, ax)

        if hist:
            npoints = self._2d_hist_fc_mw(fig, ax, nbins, wave_type)
        else:
            npoints = self._scatter_fc_mw(
                fig, ax, wave_type, colorby, colormap)
        if fit:
            self._fit_fc_mw(vel, k, ax, slope=slope)
        extra_text = 'Stress drop curves'
        self._set_plot_title(ax, npoints, extra_text)
        self._add_grid(ax_Mo)
        ax_Mo.set_ylabel(get_param_label('fc'))
        plt.show()

    def plot_Er_mw(self, hist=False, fit=False, nbins=None,
                   colorby=None, colormap=None):
        """
        Plot the logarithm of the radiated energy vs the moment magnitude.

        Parameters
        ----------
        hist : bool
            If True, plot a 2D histogram instead of a scatter plot.
        fit : bool
            If True, plot a linear regression of Er vs mw.
        nbins : int
            Number of bins for the 2D histogram.
        colorby : str
            Color the data points by the given parameter.
        colormap : str
            Use this colormap for the colorby parameter instead of the default.
        """
        fig, ax, ax_Mo = self._make_mw_axis()
        ax_Mo.set_yscale('log')

        mu = np.nanmean((self.vs * 1e3)**2 * self.rho)
        self._apparent_stress_curves_Er_mw(mu, ax)

        if hist:
            npoints = self._2d_hist_Er_mw(fig, ax, nbins)
        else:
            npoints = self._scatter_Er_mw(fig, ax, colorby, colormap)
        if fit:
            raise NotImplementedError('Fit not implemented yet for Er_mw')

        extra_text = 'Apparent stress curves'
        self._set_plot_title(ax, npoints, extra_text)
        self._add_grid(ax_Mo)
        ax_Mo.set_ylabel(get_param_label('Er'))
        plt.show()

    def plot_ssd_mw(self, hist=False, fit=False, nbins=None,
                    colorby=None, colormap=None):
        """
        Plot the logarithm of static stress drop vs moment magnitude.

        Parameters
        ----------
        hist : bool
            If True, plot a 2D histogram instead of a scatter plot.
        fit : bool
            If True, plot a linear regression of ssd vs mw.
        nbins : int
            Number of bins for the 2D histogram.
        colorby : str
            Color the data points by the given parameter.
        colormap : str
            Use this colormap for the colorby parameter instead of the default.
        """
        fig, ax, ax_Mo = self._make_mw_axis()
        ax_Mo.set_yscale('log')

        ssd_min = np.min(self.ssd - self.ssd_err_minus)
        ssd_max = np.max(self.ssd + self.ssd_err_plus)
        ssd_min = 10**(np.floor(np.log10(ssd_min)))
        ssd_max = 10**(np.ceil(np.log10(ssd_max)))
        ax.set_ylim(ssd_min, ssd_max)

        if hist:
            self._2d_hist_ssd_mw(fig, ax, nbins)
        else:
            self._scatter_ssd_mw(fig, ax, colorby, colormap)
        if fit:
            raise NotImplementedError('Fit not implemented yet for ssd_mw')

        self._set_plot_title(ax)
        self._add_grid(ax_Mo)
        ax_Mo.set_ylabel(get_param_label('ssd'))
        plt.show()

    def plot_ssd_depth(self, hist=False, fit=False, nbins=None,
                       colorby=None, colormap=None):
        """
        Plot the logarithm of static stress drop vs depth.

        Parameters
        ----------
        hist : bool
            If True, plot a 2D histogram instead of a scatter plot.
        fit : bool
            If True, plot a linear regression of ssd vs depth.
        nbins : int
            Number of bins for the 2D histogram.
        colorby : str
            Color the data points by the given parameter.
        colormap : str
            Use this colormap for the colorby parameter instead of the default.
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        depth_min = np.min(self.depth)
        depth_max = np.max(self.depth)
        padding = 0.15 * (depth_max - depth_min)
        ax.set_xlim(depth_min-padding, depth_max+padding)
        ax.set_xlabel(get_param_label('depth'))
        ax.set_ylabel(get_param_label('ssd'))
        ax.set_yscale('log')
        ax.set_ylim(1e-3, 1e3)

        if hist:
            self._2d_hist_ssd_depth(fig, ax, nbins)
        else:
            self._scatter_ssd_depth(fig, ax, colorby, colormap)
        if fit:
            raise NotImplementedError(
                'Fit not implemented yet for ssd_depth')

        self._set_plot_title(ax)
        self._add_grid(ax)
        plt.show()

    def plot_mw_time(self, colorby=None, colormap=None):
        """
        Plot the moment magnitude vs time.

        Parameters
        ----------
        colorby : str
            Color the data points by the given parameter.
        colormap : str
            Use this colormap for the colorby parameter instead of the default.
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.set_zorder(10)
        ax.patch.set_visible(False)
        ax.set_xlabel('Time')
        ax.xaxis.set_major_formatter(
            plt.matplotlib.dates.DateFormatter('%Y-%m-%d')
        )
        ax.set_ylabel(get_param_label('Mw'))
        mpl_time = np.array([t.matplotlib_date for t in self.time])
        padding_x = 0.05 * (np.max(mpl_time) - np.min(mpl_time))
        padding_y = 0.05 * (
            np.max(self.Mw + self.Mw_err_plus) -
            np.min(self.Mw - self.Mw_err_minus)
        )
        ax.set_xlim(np.min(mpl_time) - padding_x, np.max(mpl_time) + padding_x)
        ax.set_ylim(
            np.min(self.Mw - self.Mw_err_minus) - padding_y,
            np.max(self.Mw + self.Mw_err_plus) + padding_y
        )
        cond = ~np.isnan(self.Mw)
        npoints = self._scatter_plot(
            fig, ax, mpl_time, self.Mw,
            np.zeros_like(self.time), [self.Mw_err_minus, self.Mw_err_plus],
            self.evids, 'Mw {:.2f}', cond,
            colorby, colormap, cbaxes_location='right', xformat='{}')

        # Add cumulative moment axis
        ax2 = ax.twinx()
        ax2.set_zorder(0)
        ax2.set_ylabel('Cumulative Moment (N·m)')
        isort = np.argsort(mpl_time)
        mpl_time = mpl_time[isort]
        Mo = self.Mo[isort]
        cumulative_moment = np.cumsum(Mo)
        ax2.plot(
            mpl_time, cumulative_moment, color='red', alpha=0.6, lw=2)
        ax2.set_ylim(0, np.max(cumulative_moment) * 1.1)
        final_moment = cumulative_moment[-1]
        cum_mag = moment_to_mag(final_moment)
        cum_mag_text = f'Cumulative\nMw: {cum_mag:.2f}'
        bbox = {
            'facecolor': 'white',
            'edgecolor': 'black',
            'boxstyle': 'round',
            'alpha': 0.9
        }
        # we add the text to the first axis, so that it is on top,
        # but we need to transform the coordinates to the second axis
        ax.text(
            mpl_time[-1], cumulative_moment[-1],
            cum_mag_text, color='red',
            transform=ax2.transData,
            verticalalignment='bottom', horizontalalignment='left',
            bbox=bbox, zorder=50)

        extra_text = 'Mw vs Time'
        self._set_plot_title(ax, npoints, extra_text, align='left')
        self._add_grid(ax)
        plt.show()

    def plot_hist(self, param_name, nbins=None, wave_type='S'):
        """
        Plot a histogram of the given parameter.

        Parameters
        ----------
        param_name : str
            Name of the parameter to plot.
        nbins : int
            Number of bins for the histogram.
        wave_type : str
            One of "P", "S", "SV" or "SH". Default is "S".
        """
        _, unit, log, color = valid_parameters[param_name]
        log = log == 'log'
        values = getattr(self, param_name)
        if param_name == 'fc':
            values = values[self.wave_type == wave_type]
            if len(values) == 0:
                raise ValueError(
                    f'No events found for wave type "{wave_type}"')
        values = values[~np.isnan(values)]
        nvalues = len(values)
        if nvalues < 2:
            raise ValueError(f'Not enough values to plot histogram: {nvalues}')
        if nbins is None:
            # find the number of bins as the closest power of 10 to the
            # number of values divided by 10
            nbins_exp = int(np.ceil(np.log10(nvalues / 10)))
            nbins = int(10**nbins_exp)
            # reduce the number of bins if there are too few values per bin
            if nvalues / nbins < 2:
                nbins //= 2
        if log:
            values_mean = 10**(np.mean(np.log10(values)))
            min_exp = np.floor(np.min(np.log10(values)))
            max_exp = np.ceil(np.max(np.log10(values)))
            bins = np.logspace(min_exp, max_exp, nbins + 1)
        else:
            values_mean = np.mean(values)
            maxval = np.max(np.abs(values))
            bins = np.linspace(0, maxval, nbins)
        _fig, ax = plt.subplots()
        ax.hist(values, bins=bins, color=color)
        ax.axvline(values_mean, color='red')
        txt = (
            f'  mean: {values_mean:.1e}'
            if log else
            f'  mean: {values_mean:.4f}'
        )
        if unit is not None:
            txt += f' {unit}'
        ax.text(
            values_mean, 0.95, txt, color='red',
            transform=ax.get_xaxis_transform())
        if log:
            ax.set_xscale('log')
        ax.set_xlabel(get_param_label(param_name))
        ax.set_ylabel('Counts')
        title = f'{nvalues} values, {nbins} bins - {self.stat}'
        if self.runid is not None:
            title += f' - runid {self.runid}'
        ax.set_title(title)
        plt.show()

    def plot_gutenberg_richter(
            self, magmin=None, magmax=None, nbins=None, fit=False):
        """
        Plot the Gutenberg-Richter law.

        Parameters
        ----------
        magmin : float
            Minimum magnitude for the plot.
        magmax : float
            Maximum magnitude for the plot.
        nbins : int
            Number of bins for the histogram.
        fit : bool
            If True, find the magnitude of completeness and
            fit the Gutenberg-Richter law to the data.
        """
        if magmin is None:
            magmin = np.floor(np.min(self.Mw - self.Mw_err_minus))
        if magmax is None:
            magmax = np.ceil(np.max(self.Mw + self.Mw_err_plus))
        if nbins is None:
            nbins = int((magmax - magmin) * 10)
            print(f'Number of bins autoset to: {nbins}')

        bins = np.linspace(magmin, magmax, nbins + 1)
        hist, bin_edges = np.histogram(self.Mw, bins=bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        cum_nevs = np.cumsum(hist[::-1])[::-1]

        if fit:
            mc = compute_mc(self.Mw, 0.1, bin_centers, 0.1)
            print(f'Magnitude of completeness: {mc:.2f}')
            a, b = fit_gutenberg_richter(bin_centers, cum_nevs, mc)
            print(f'G-R law parameters: a={a:.2f}, b={b:.2f}')
            mw_fit = np.linspace(mc, magmax, 100)
            gr_fit = 10**(a - b * mw_fit)

        _fig, ax1 = plt.subplots(figsize=(10, 6))
        ax1.set_zorder(1)
        ax1.patch.set_visible(False)
        ax1.plot(
            bin_centers, cum_nevs, 'o',
            color='#FCBA25', markeredgecolor='black', markersize=8,
            zorder=10)
        if fit:
            outline = patheffects.withStroke(linewidth=0.5, foreground='black')
            bbox = {
                'facecolor': 'white',
                'edgecolor': 'black',
                'boxstyle': 'round',
                'alpha': 0.7
            }
            # plot the G-R law
            ax1.plot(mw_fit, gr_fit, color='#FCBA25', linewidth=4, zorder=2)
            xmin, xmax = ax1.get_xlim()
            # add a label with the G-R law parameters
            gr_text_x = 0.805
            gr_text_x_data = xmin + gr_text_x * (xmax - xmin)
            gr_text_y_data = 10**(a - b * gr_text_x_data)*1.2
            gr_text = f'G-R params:\na={a:.2f}, b={b:.2f}'
            ax1.text(
                gr_text_x, gr_text_y_data, gr_text, color='#FCBA25',
                transform=ax1.get_yaxis_transform(),
                horizontalalignment='left', verticalalignment='bottom',
                path_effects=[outline,], bbox=bbox
            )
            # add a vertical line at the magnitude of completeness
            ax1.axvline(mc, color='red', linestyle='--', zorder=2)
            # add a label with the magnitude of completeness
            mc_text_x_data = mc + 0.02 * (xmax - xmin)
            mc_text = f'$M_C$={mc:.2f}'
            ax1.text(
                mc_text_x_data, 0.97, mc_text, color='red',
                transform=ax1.get_xaxis_transform(),
                horizontalalignment='left', verticalalignment='top',
                path_effects=[outline,], bbox=bbox
            )
        ax1.set_yscale('log')
        ax1.set_xlabel('Magnitude (Mw)')
        ax1.set_ylabel('Cumulative Number of Events', color='#FCBA25')
        ax1.tick_params(axis='y', labelcolor='#FCBA25')

        ax2 = ax1.twinx()
        ax2.set_zorder(0)
        ax2.hist(self.Mw, bins=bins, color='#1E90FF', alpha=0.6)
        ax2.set_ylabel('Number of Events', color='#1E90FF')
        ax2.tick_params(axis='y', labelcolor='#1E90FF')
        ax2.set_ylim(0, max(hist) * 1.5)

        extra_text = f'Gutenberg-Richter Law\n{nbins} bins'
        self._set_plot_title(ax1, extra_text=extra_text, align='right')
        self._add_grid(ax1)
        plt.show()


def run():
    """Run the script."""
    args = parse_args()
    params = Params(args)

    if os.path.exists('problems.txt'):
        _skip_events(params)
    params.filter(
        stamin=args.stamin, stamax=args.stamax,
        magmin=args.magmin, magmax=args.magmax,
        ssdmin=args.ssdmin, ssdmax=args.ssdmax,
        sigmaamin=args.sigmaamin, sigmaamax=args.sigmaamax,
    )
    if len(params.evids) == 0:
        raise ValueError('No events found')

    if args.plot_type == 'fc_mw':
        params.plot_fc_mw(
            args.hist, args.fit, args.slope, args.nbins, args.wave_type,
            args.colorby, args.colormap)
    elif args.plot_type == 'Er_mw':
        params.plot_Er_mw(
            args.hist, args.fit, args.nbins,
            args.colorby, args.colormap)
    elif args.plot_type == 'ssd_mw':
        params.plot_ssd_mw(
            args.hist, args.fit, args.nbins,
            args.colorby, args.colormap)
    elif args.plot_type == 'ssd_depth':
        params.plot_ssd_depth(
            args.hist, args.fit, args.nbins,
            args.colorby, args.colormap)
    elif args.plot_type == 'mw_time':
        params.plot_mw_time(args.colorby, args.colormap)
    elif args.plot_type == 'GR':
        params.plot_gutenberg_richter(
            args.magmin, args.magmax, args.nbins, args.fit)
    elif args.plot_type in valid_plot_types:
        params.plot_hist(args.plot_type, args.nbins, args.wave_type)


def _skip_events(params):
    evid_skip = np.loadtxt('problems.txt', usecols=(0,), dtype=str)
    # case in which there is only one evid:
    if not evid_skip.shape:
        evid_skip = evid_skip[None]
    skip_idx = []
    for evid in evid_skip:
        with contextlib.suppress(Exception):
            skip_idx.append(np.where(params.evids == evid)[0][0])
    skipped_evid_str = ' '.join(params.evids[skip_idx])
    print(f'Skipping events: {skipped_evid_str}')
    params.skip_events(skip_idx)


def main():
    """Main function."""
    try:
        run()
    except KeyboardInterrupt:
        sys.exit(0)
    except Exception as msg:
        sys.exit(msg)


if __name__ == '__main__':
    main()
