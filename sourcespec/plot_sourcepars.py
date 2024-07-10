# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
plot_sourcepars.py

1D or 2D plot of source parameters from a sqlite parameter file.

:copyright:
    2023-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import argparse
import contextlib
import sys
import os
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

valid_plot_types = [
    'fc', 'Er', 'ssd', 'ra', 'Mo', 't_star', 'Qo', 'sigma_a',
    'fc_mw', 'Er_mw', 'ssd_mw']


class Annot():
    """
    Annotate the plot with the evid, Mw and the value of the parameter.
    """
    def __init__(self, xdata, ydata, labels, yformat):
        self.xdata = xdata
        self.ydata = ydata
        self.labels = labels
        self.yformat = yformat

    def __call__(self, event):
        ind = event.ind
        x = np.take(self.xdata, ind)
        y = np.take(self.ydata, ind)
        label = np.take(self.labels, ind)
        for Mw, yvalue, evid in zip(x, y, label):
            ystring = self.yformat.format(yvalue)
            print(f'{evid} Mw {Mw:.1f} {ystring}')


def mag_to_moment(mag, b=0.5):
    """
    Convert magnitude to moment.

    The parameter b is used to change the slope of the fc-Mw stress drop curve.

    The standard value of b is 0.5, which corresponds to a self-similar model.
    """
    mag = np.asarray(mag)
    return np.power(10, (3 * b * mag + 9.1))


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
        self.mw = query_event_params_into_numpy(
            self.cur, f'Mw_{stat}', np.float64, query_condition)
        self.mw_err_minus = query_event_params_into_numpy(
            self.cur, f'Mw_{stat}_err_minus', np.float64, query_condition)
        self.mw_err_plus = query_event_params_into_numpy(
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
        self.vp = np.delete(self.vp, idx)
        self.vs = np.delete(self.vs, idx)
        self.rho = np.delete(self.rho, idx)
        self.kp = np.delete(self.kp, idx)
        self.ks = np.delete(self.ks, idx)
        self.wave_type = np.delete(self.wave_type, idx)
        self.nsta = np.delete(self.nsta, idx)
        self.Mo = np.delete(self.Mo, idx)
        self.mw = np.delete(self.mw, idx)
        self.mw_err_minus = np.delete(self.mw_err_minus, idx)
        self.mw_err_plus = np.delete(self.mw_err_plus, idx)
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
            cond = np.logical_and(cond, self.mw >= magmin)
        if magmax is not None:
            cond = np.logical_and(cond, self.mw <= magmax)
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
            mag_min = np.floor(np.min(self.mw - self.mw_err_minus))
        if self.args.magmax is not None:
            mag_max = self.args.magmax
        else:
            mag_max = np.ceil(np.max(self.mw + self.mw_err_plus))
        xlim_mag = (mag_min, mag_max)
        fig = plt.figure(figsize=(10, 6))
        ax_Mo = fig.add_subplot(111)
        ax_Mo.set_xscale('log')
        ax_Mo.set_xlim(mag_to_moment(xlim_mag))
        ax = ax_Mo.twiny()
        ax.set_xlim(xlim_mag)
        ax.set_xlabel('Mw')
        ax_Mo.set_xlabel('Mo (N·m)')
        return fig, ax, ax_Mo

    def _add_grid(self, ax):
        """Add grid to the plot."""
        ax.grid(
            True, which='major', linestyle='solid',
            color='#999999', zorder=0)
        ax.grid(
            True, which='minor', linestyle='solid',
            color='#DDDDDD', zorder=0)

    def _set_plot_title(self, ax, nevs=None, extra_text=None):
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
        ax.set_title(title, y=0.95, verticalalignment='top')

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
            ax.text(mw_test[-1], fc_test[-1], label)
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
            ax.text(mw_test[-1], Er_test[-1], label)
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
            ax.text(mw_test[-1], Er_test[-1], label)
        ax.set_ylim((Er_min * 0.5, Er_max * 2))

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
        mw = self.mw[cond]
        if len(mw) == 0:
            raise ValueError(f'No events found for wave type "{wave_type}"')
        fc = self.fc[cond]
        counts, _, _ = np.histogram2d(mw, fc, bins=(mw_bins, fc_bins))
        cm = ax.pcolormesh(
            mw_bins[:-1], fc_bins[:-1], counts.T,
            cmap='magma_r', shading='auto')
        cbaxes = fig.add_axes([0.15, 0.15, 0.02, 0.2])
        plt.colorbar(cm, cax=cbaxes, orientation='vertical', label='counts')
        return len(fc)

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
        # Er can be NaN
        cond = ~np.isnan(self.Er)
        mw = self.mw[cond]
        Er = self.Er[cond]
        counts, _, _ = np.histogram2d(mw, Er, bins=(mw_bins, Er_bins))
        cm = ax.pcolormesh(
            mw_bins[:-1], Er_bins[:-1], counts.T,
            cmap='magma_r', shading='auto')
        cbaxes = fig.add_axes([0.15, 0.15, 0.02, 0.2])
        plt.colorbar(cm, cax=cbaxes, orientation='vertical', label='counts')
        return len(Er)

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
        counts, _, _ = np.histogram2d(
            self.mw, self.ssd, bins=(mw_bins, ssd_bins))
        cm = ax.pcolormesh(
            mw_bins[:-1], ssd_bins[:-1], counts.T,
            cmap='magma_r', shading='auto')
        cbaxes = fig.add_axes([0.15, 0.15, 0.02, 0.2])
        plt.colorbar(cm, cax=cbaxes, orientation='vertical', label='counts')

    def _scatter_fc_mw(self, fig, ax, wave_type='S'):
        """Plot the scatter plot of fc vs mw."""
        cond = self.wave_type == wave_type
        evids = self.evids[cond]
        if len(evids) == 0:
            raise ValueError(f'No events found for wave type "{wave_type}"')
        mw = self.mw[cond]
        mw_err_minus = self.mw_err_minus[cond]
        mw_err_plus = self.mw_err_plus[cond]
        fc = self.fc[cond]
        fc_err_minus = self.fc_err_minus[cond]
        fc_err_plus = self.fc_err_plus[cond]
        alpha = 1
        ax.errorbar(
            mw, fc,
            xerr=[mw_err_minus, mw_err_plus],
            yerr=[fc_err_minus, fc_err_plus],
            fmt='o', mec='black', mfc='#FCBA25', ecolor='#FCBA25',
            alpha=alpha)
        ax.scatter(mw, fc, alpha=0, picker=True, zorder=20)
        yformat = 'fc {:.2f} Hz'
        annot = Annot(mw, fc, evids, yformat)
        fig.canvas.mpl_connect('pick_event', annot)
        return len(fc)

    def _scatter_Er_mw(self, fig, ax):
        """Plot the scatter plot of Er vs mw."""
        # Er can be NaN
        cond = ~np.isnan(self.Er)
        mw = self.mw[cond]
        mw_err_minus = self.mw_err_minus[cond]
        mw_err_plus = self.mw_err_plus[cond]
        Er = self.Er[cond]
        Er_err_minus = self.Er_err_minus[cond]
        Er_err_plus = self.Er_err_plus[cond]
        evids = self.evids[cond]
        alpha = 1
        ax.errorbar(
            mw, Er,
            xerr=[mw_err_minus, mw_err_plus],
            yerr=[Er_err_minus, Er_err_plus],
            fmt='o', mec='black', mfc='#FCBA25', ecolor='#FCBA25',
            alpha=alpha)
        ax.scatter(mw, Er, alpha=0, picker=True, zorder=20)
        yformat = 'Er {:.1e} N·m'
        annot = Annot(mw, Er, evids, yformat)
        fig.canvas.mpl_connect('pick_event', annot)
        return len(Er)

    def _scatter_ssd_mw(self, fig, ax):
        """Plot the scatter plot of ssd vs mw."""
        alpha = 1
        ax.errorbar(
            self.mw, self.ssd,
            xerr=[self.mw_err_minus, self.mw_err_plus],
            yerr=[self.ssd_err_minus, self.ssd_err_plus],
            fmt='o', mec='black', mfc='#FCBA25', ecolor='#FCBA25',
            alpha=alpha)
        ax.scatter(self.mw, self.ssd, alpha=0, picker=True, zorder=20)
        yformat = 'ssd {:.2e} MPa'
        annot = Annot(self.mw, self.ssd, self.evids, yformat)
        fig.canvas.mpl_connect('pick_event', annot)

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
        popt, _pcov = curve_fit(f, self.mw, y)
        try:
            a, b = popt
        except Exception:
            a, = popt
            b = -0.5
        print(f'a: {a:.1f} b {b:.1f}:')
        slope = (3 / 2) / b
        print(f'slope: {slope:.1f}')
        r2 = calc_r2(self.mw, y, np.ones_like(y), a, b)
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
                   wave_type='S'):
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
            npoints = self._scatter_fc_mw(fig, ax, wave_type)
        if fit:
            self._fit_fc_mw(vel, k, ax, slope=slope)
        extra_text = 'Stress drop curves'
        self._set_plot_title(ax, npoints, extra_text)
        self._add_grid(ax_Mo)
        ax_Mo.set_ylabel('fc (Hz)')
        plt.show()

    def plot_Er_mw(self, hist=False, fit=False, nbins=None):
        """
        Plot the logarithm of the radiated energy vs the moment magnitude.

        Parameters
        ----------
        hist : bool
            If True, plot a 2D histogram instead of a scatter plot.
        fit : bool
            If True, plot a linear regression of Er vs mw.
        slope : bool
            If True, also fit the slope of the linear regression.
        """
        fig, ax, ax_Mo = self._make_mw_axis()
        ax_Mo.set_yscale('log')

        mu = np.nanmean((self.vs * 1e3)**2 * self.rho)
        self._apparent_stress_curves_Er_mw(mu, ax)

        if hist:
            npoints = self._2d_hist_Er_mw(fig, ax, nbins)
        else:
            npoints = self._scatter_Er_mw(fig, ax)
        if fit:
            raise NotImplementedError('Fit not implemented yet for Er_mw')

        extra_text = 'Apparent stress curves'
        self._set_plot_title(ax, npoints, extra_text)
        self._add_grid(ax_Mo)
        ax_Mo.set_ylabel('Er (N·m)')
        plt.show()

    def plot_ssd_mw(self, hist=False, fit=False, nbins=None):
        """
        Plot the logarithm of static stress drop vs moment magnitude.

        Parameters
        ----------
        hist : bool
            If True, plot a 2D histogram instead of a scatter plot.
        fit : bool
            If True, plot a linear regression of ssd vs mw.
        slope : bool
            If True, also fit the slope of the linear regression.
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
            self._scatter_ssd_mw(fig, ax)
        if fit:
            raise NotImplementedError('Fit not implemented yet for ssd_mw')

        self._set_plot_title(ax)
        self._add_grid(ax_Mo)
        ax_Mo.set_ylabel('ssd (MPa)')
        plt.show()

    def plot_hist(self, param_name, nbins=None, wave_type='S'):
        """Plot a histogram of the given parameter."""
        parameters = {
            'fc': ('Corner Frequency', 'Hz', 'log'),
            'ssd': ('Static Stress Drop', 'MPa', 'log'),
            'ra': ('Source Radius', 'm', 'lin'),
            'Mo': ('Seismic Moment', 'N·m', 'log'),
            't_star': ('T star', 's', 'lin'),
            'Qo': ('Qo', None, 'lin'),
            'Er': ('Radiated Energy', 'N·m', 'log'),
            'sigma_a': ('Apparent Stress', 'MPa', 'log'),
        }
        description, unit, log = parameters[param_name]
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
        ax.hist(values, bins=bins)
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
        if unit is not None:
            ax.set_xlabel(f'{description} ({unit})')
        else:
            ax.set_xlabel(f'{description}')
        ax.set_ylabel('Counts')
        title = f'{nvalues} values, {nbins} bins - {self.stat}'
        if self.runid is not None:
            title += f' - runid {self.runid}'
        ax.set_title(title)
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
            args.hist, args.fit, args.slope, args.nbins, args.wave_type)
    elif args.plot_type == 'Er_mw':
        params.plot_Er_mw(args.hist, args.fit, args.nbins)
    elif args.plot_type == 'ssd_mw':
        params.plot_ssd_mw(args.hist, args.fit, args.nbins)
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
