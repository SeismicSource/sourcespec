# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Object-oriented python implementation of the Fortran code in the
TRVDRV subroutine of hypo71 to calculate travel times and takeoff
and incidence angles for P and S waves using a 1D velocity model.

:copyright:
    2020-2025 Kris Vanneste <kris.vanneste@oma.be>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import numpy as np
from scipy.optimize import minimize_scalar


class CrustalVelocityModel:
    """
    1D layered velocity model

    :param depths:
        depth of layer tops (in km)
    :param VP:
        P-wave veolicities (in km/s)
    :param VS:
        S-wave velocities (in km/s)
    :param name:
        str, model name
    """
    def __init__(self, depths, VP, VS, name=''):
        assert len(depths) == len(VP) == len(VS)
        self.depths = np.asarray(depths, dtype='f')
        self.VP = np.asarray(VP, dtype='f')
        self.VS = np.asarray(VS, dtype='f')
        self.name = name
        # Layer thicknesses (excluding halfspace)
        self.thicknesses = np.diff(depths)
        # Refraction matrices
        self.RXTT = None
        self.RXHD = None

    def __len__(self):
        return len(self.depths)

    def __repr__(self):
        return f'<CrustalVelocityModel "{self.name}" (n={self.num_layers})>'

    @property
    def num_layers(self):
        """
        Number of layers, excluding halfspace
        """
        return len(self) - 1

    @property
    def max_depth(self):
        """
        Maximum depth
        """
        return self.depths[-1]

    def recalc_vs_from_vp(self, vp_vs_ratio):
        """
        Compute and set VS based on VP and VP/VS ratio

        :param vp_vs_ratio:
            float, VP/VS ratio
        """
        self.VS = self.VP / vp_vs_ratio

    def get_velocities(self, wave='P'):
        """
        Get velocities corresponding to given wave type

        :param wave:
            char ('P' or 'S') or int (0 = P, 1 = S)
            (default: 'P')

        :return:
            1D array, velocities in km/s [NL + 1]
        """
        if wave == 0:
            V = self.VP
        elif wave == 1:
            V = self.VS
        elif wave.upper() == 'P':
            V = self.VP
        elif wave.upper() == 'S':
            V = self.VS
        else:
            raise ValueError(f'Unknown wave type: {wave}')
        return V

    def get_velocity_ratios(self, wave='P'):
        """
        Compute ratio between layer velocity and velocity in layer below

        :param wave:
            see :meth:`get_velocities`

        :return:
            1D array [NL]
        """
        V = self.get_velocities(wave=wave)
        return V[:-1] / V[1:]

    def get_velocity_at_depth(self, depth, wave='P'):
        """
        Get velocity at particular depth

        :param depth:
            float, depth (in km)
        :param wave:
            see :meth:`get_velocities`

        :return:
            float, velocity (in km/s)
        """
        idx = self.get_layer_index(depth)
        return self.get_velocities(wave=wave)[idx]

    def calc_critical_angles(self, wave='P', deg=False):
        """
        Compute critical angles in each layer for given wave type

        :param wave:
            see :meth:`get_velocities`
        :param deg:
            bool, whether or not to return angles as degrees
            (default: False, will return as radians)

        :return:
            1D array, critical angles [NL]
        """
        angles = np.arcsin(self.get_velocity_ratios(wave=wave))
        if deg:
            angles = np.degrees(angles)
        return angles

    def calc_travel_times(self, layer_angles, wave='P'):
        """
        Compute travel times in each layer corresponding to given
        layer angles and wave type

        :param layer_angles:
            float or 1D array, angles in radians [NL or less]
            If array and length is smaller than number of layers in model,
            bottom layers are supposed to be missing, and the result
            will be truncated to the same length
        :param wave:
            see :meth:`get_velocities`

        :return:
            1D array, travel times in seconds, with length:
            - NL if :param:`layer_angles` is scalar
            - length of layer_angles array.
        """
        NL = (
            self.num_layers if np.isscalar(layer_angles)
            else len(layer_angles)
        )
        V = self.get_velocities(wave=wave)[:NL]
        THK = self.thicknesses[:NL]
        DIS = THK / np.cos(layer_angles)
        return DIS / V

    def calc_horizontal_distances(self, layer_angles):
        """
        Compute horizontal distances in each layer corresponding to
        given layer angles

        :param layer_angles:
            see :meth:`calc_travel_times`

        :return:
            1D array, horizontal distances in km, with length:
            - NL if :param:`layer_angles` is scalar
            - length of layer_angles array.
        """
        NL = len(layer_angles)
        THK = self.thicknesses[:NL]
        return THK * np.tan(layer_angles)

    def calc_layer_angles_above(
            self, bottom_layer_idx, bottom_angle, wave='P'):
        """
        Compute layer angles in higher layers for a ray with a given
        angle in the layer below

        :param bottom_layer_idx:
            int, index of bottom layer
        :param bottom_angle:
            float, angle in bottom layer (in radians)
        :param wave:
            see :meth:`get_velocities`

        :return:
            1D array, with length equal to :param:`bottom_layer_idx` + 1,
            i.e., bottom layer is included
        """
        if bottom_layer_idx < 0:
            bottom_layer_idx = self.num_layers + bottom_layer_idx
        assert bottom_layer_idx <= self.num_layers

        # TODO: should we check if angles remain below critical angle
        # in each layer?
        Vratios = self.get_velocity_ratios(wave=wave)[:bottom_layer_idx]
        layer_angles = [bottom_angle]
        for L in range(bottom_layer_idx - 1, -1, -1):
            angle = np.arcsin(Vratios[L] * np.sin(layer_angles[-1]))
            layer_angles.append(angle)

        return np.array(layer_angles[::-1])

    def get_layer_index(self, Z):
        """
        Find index of layer containing a particular depth

        Note: if Z is smaller than top depth in model, returned index
        is zero.

        :param Z:
            float, depth (in km)

        :return:
            int, layer index
        """
        try:
            return np.where(self.depths <= Z)[0][-1]
        except IndexError:
            # Z is above top depth (e.g., station at altitude > 0), return 0
            return 0

    def calc_refraction_matrices(self, recalc=False):
        """
        Compute cumulative travel times and horizontal distances
        for downgoing rays from any top layer critically refracting
        in any bottom layer

        :param recalc:
            bool, whether or not to force recalculating the matrices
            (default: False)

        :return:
            (RXTT, RXHD) tuples of 3D arrays [2, NL, NL], with dimensions
            corresponding to [wave types, top layers, bottom layers]
            - RXTT: cumulative travel times
            - RXHD: cumulative horizontal distances
            Note: these matrices are calculated only once, and stored in
            :prop:`RXTT` and :prop:`RXHD`
        """
        if self.RXTT is None or self.RXHD is None or recalc:
            NL = self.num_layers
            self.RXTT = np.zeros((2, NL, NL))
            self.RXHD = np.zeros((2, NL, NL))
            for W, wave in enumerate(('P', 'S')):
                theta_crit = self.calc_critical_angles(wave=wave)
                for BL in range(NL):
                    bottom_angle = theta_crit[BL]
                    layer_angles = self.calc_layer_angles_above(
                        BL, bottom_angle, wave=wave)
                    travel_times = self.calc_travel_times(
                        layer_angles, wave=wave)
                    cumul_travel_times = np.cumsum(travel_times[::-1])[::-1]
                    hdistances = self.calc_horizontal_distances(layer_angles)
                    cumul_hdistances = np.cumsum(hdistances[::-1])[::-1]
                    self.RXTT[W, :BL+1, BL] = cumul_travel_times
                    self.RXHD[W, :BL+1, BL] = cumul_hdistances
        return (self.RXTT, self.RXHD)

    def calc_downgoing_ray_tt_and_hd(self, Z, wave='P'):
        """
        Compute cumulative travel times / horizontal distances for
        a downgoing ray originating at depth Z, critically refracting
        at different bottom layers.

        :param Z:
            float, depth at which ray originates (in km)
        :param wave:
            see :meth:`get_velocities`

        :return:
            (rxtt, rxhd) tuple of 1D arrays [NL]
            - rxtt: cumulative travel times (in seconds)
            - rxhd: cumulative horizontal distances (in km)
            Note: values in layers above Z are zero
        """
        RXTT, RXHD = self.calc_refraction_matrices()
        if wave in (0, 1):
            W = wave
        elif wave.upper() == 'P':
            W = 0
        elif wave.upper() == 'S':
            W = 1
        else:
            raise ValueError(f'Unknown wave type: {wave}')
        RXTT, RXHD = RXTT[W], RXHD[W]
        TL = self.get_layer_index(Z)
        NL = self.num_layers
        if TL >= NL:
            rxtt = np.zeros(NL) * np.nan
            rxhd = np.zeros(NL) * np.nan
        else:
            # Correct for depth inside layer
            # (works also for negative depths inside layer, which are treated
            # as lying above the layer, but with the same velocities)
            vdl = Z - self.depths[TL]
            thkl = self.thicknesses[TL]
            # Note: take copy, otherwise original RXTT and RXHD arrays
            # are modified!
            rxtt, rxhd = RXTT[TL].copy(), RXHD[TL].copy()
            if TL < NL-1:
                # Time / distance in layer in which Z is situated
                # (for different bottom layers)
                ttl = RXTT[TL] - RXTT[TL+1]
                hdl = RXHD[TL] - RXHD[TL+1]
            elif TL == NL - 1:
                ttl, hdl = RXTT[TL], RXHD[TL]
            else:
                raise RuntimeError(
                    f'Invalid layer index: TL={TL} >= NL={NL}. '
                    f'Z={Z} exceeds model depth range [0, {self.max_depth}].'
                )
            idxs = rxtt > 0
            rxtt[idxs] -= (ttl[idxs] * vdl/thkl)
            rxhd[idxs] -= (hdl[idxs] * vdl/thkl)
        return (rxtt, rxhd)

    def calc_refwav_tt(self, Zf, Zs, Repi, wave='P'):
        """
        Compute travel times for refracted waves in different bottom
        layers between focus at depth Zf and station at depth Zs,
        separated by epicentral distance Repi.

        :param Zf:
            float, depth of focus (in km)
        :param Zs:
            float, depth of station (in km)
        :param Repi:
            float, epicentral distance (= horizontal distance between
            earthquake focus and receiver station
        :param wave:
            see :meth:`get_velocities`

        :return:
            1D array, travel times corresponding to different bottom
            layers (in seconds) [NL]
            Note: travel times in layers above the maximum of Zf and Zs
            are set to NaN. If both Zf and Zs are situated in the half-
            space, only NaN values will be returned.
        """
        rxttf, rxhdf = self.calc_downgoing_ray_tt_and_hd(Zf, wave=wave)
        rxtts, rxhds = self.calc_downgoing_ray_tt_and_hd(Zs, wave=wave)
        V = self.get_velocities(wave)

        NL = self.num_layers
        Lf = self.get_layer_index(Zf)
        Ls = self.get_layer_index(Zs)
        # highest possible bottom layer
        BLmin = max(Lf, Ls)

        # Limit BLmin to NL - 1
        # if np.isclose(max(Zf, Zs), self.max_depth):
        #    BLmin = NL - 1
        # Capture case where lower of Zf or Zs corresponds to layer interface
        if np.isclose(max(Zf, Zs), self.depths).any():
            BLmin -= 1

        refwav_tt = np.zeros(NL) * np.nan
        if BLmin < NL:
            # If Zf and/or Zs are in halfspace, there is no refracted wave
            # if not (Lf == NL and Ls == NL):
            for BL in range(BLmin, NL):
                VTT = rxttf[BL] + rxtts[BL]
                # Distance traveled horizontally along top of layer below
                HD = Repi - (rxhdf[BL] + rxhds[BL])
                if HD >= 0:
                    HTT = HD / V[BL+1]
                    refwav_tt[BL] = VTT + HTT

        return refwav_tt

    def calc_min_refwav_tt(self, Zf, Zs, Repi, wave='P'):
        """
        Compute minimum travel time for a refracted wave between focus
        at depth Zf and station at depth Zs, separated by epicentral
        distance Repi.

        :param Zf:
        :param Zs:
        :param Repi:
        :param wave:
            see :meth:`calc_refwav_tt`

        :return:
            (BL, tmin) tuple:
            - BL: int, bottom layer index
            - tmin: float, minimum travel time in seconds
            May be (None, np.inf) if refraction is not possible
        """
        refwav_tt = self.calc_refwav_tt(Zf, Zs, Repi, wave=wave)
        if not np.isnan(refwav_tt).all():
            BL = np.nanargmin(refwav_tt)
            tmin = refwav_tt[BL]
        else:
            BL, tmin = None, np.inf
        return (BL, tmin)

    def constrain_depth_range(self, Zh, Zl):
        """
        Constrain model between upper and lower depth

        Notes:
        - if Zh is smaller than top depth of first layer,
        thickness of top layer will be increased
        - if Zl is larger than top of halfspace, an additional
        layer will be added with the same velocities as the halfspace

        :param Zh:
            float, upper depth (km)
        :param Zl:
            float, lower depth (km)

        :return:
            instance of :class:`VelocityModel`
        """
        Lh = self.get_layer_index(Zh)
        Ll = self.get_layer_index(Zl)
        # Avoid bottom layer with zero thickness!
        if self.depths[Ll] == Zl:
            Ll -= 1
        # Add layer in halfspace if Zl > max_depth
        if Zl > self.max_depth:
            depths = np.hstack([self.depths, [Zl]])
            VP = np.hstack([self.VP, [self.VP[-1]]])
            VS = np.hstack([self.VS, [self.VS[-1]]])
        else:
            depths, VP, VS = self.depths, self.VP, self.VS
        depths2 = depths[Lh:Ll+2] - depths[Lh]
        depths2[1:] -= (Zh - depths[Lh])
        depths2[-1] = Zl - Zh
        VP2 = VP[Lh:Ll+2]
        VS2 = VS[Lh:Ll+2]
        return self.__class__(depths2, VP2, VS2)

    def find_reflection_ray_angles_and_tts(self, Zf, Zs, Repi, wave='P'):
        """
        Determine reflection angles and travel times of reflected rays
        in different bottom layers between focus at depth Zf and station
        at depth Zs, separated by epicentral distance Repi.

        :param Zf:
            float, depth of focus (in km)
        :param Zs:
            float, depth of station (in km)
        :param Repi:
            float, epicentral distance (= horizontal distance between
            earthquake focus and receiver station
        :param wave:
            see :meth:`get_velocities`

        Note: there is no constraint on the relative depth of Zf and Zs
        (they may be interchanged)

        :return:
            (refl_thetas, refl_tts) tuple:
            - refl_thetas: 1D array, reflection angles in radians [NL]
            - refl_tts: 1D array, travel times in seconds [NL]
        """
        # Determine higher/lower depth and layer index
        Zh = min(Zf, Zs)
        Zl = max(Zf, Zs)
        # Lh is not used, commented out to avoid linting warning
        # Lh = self.get_layer_index(Zh)
        Ll = self.get_layer_index(Zl)

        # Construct secondary velocity models from Zh and Zl downwards
        vmodel1 = self.constrain_depth_range(Zh, self.depths[-1])
        vmodel2 = self.constrain_depth_range(Zl, self.depths[-1])

        # BL = Ll - self.num_layers
        NL = self.num_layers
        refl_tts = np.zeros(NL) * np.nan
        refl_thetas = np.zeros(NL) * np.nan

        def minimize_func(theta):
            layer_angles1 = vmodel1.calc_layer_angles_above(
                BL, theta, wave=wave)
            hdistances1 = vmodel1.calc_horizontal_distances(layer_angles1)
            tot_hdistance1 = np.sum(hdistances1)
            layer_angles2 = vmodel2.calc_layer_angles_above(
                BL, theta, wave=wave)
            hdistances2 = vmodel2.calc_horizontal_distances(layer_angles2)
            tot_hdistance2 = np.sum(hdistances2)
            tot_hdistance = tot_hdistance1 + tot_hdistance2

            return np.abs(Repi - tot_hdistance)

        for BL in range(-1, Ll - NL - 1, -1):
            result = minimize_scalar(
                minimize_func, bounds=(0, np.pi/2.), method='bounded')
            if result.success:
                theta = result.x
            else:
                raise RuntimeError(result.message)
            layer_angles1 = vmodel1.calc_layer_angles_above(
                BL, theta, wave=wave)
            travel_times1 = vmodel1.calc_travel_times(layer_angles1, wave=wave)
            layer_angles2 = vmodel2.calc_layer_angles_above(
                BL, theta, wave=wave)
            travel_times2 = vmodel2.calc_travel_times(layer_angles2, wave=wave)
            tt = np.sum(travel_times1) + np.sum(travel_times2)

            refl_thetas[BL] = theta
            refl_tts[BL] = tt

        return (refl_thetas, refl_tts)

    def calc_min_reflection_tt(self, Zf, Zs, Repi, wave='P'):
        """
        Compute minimum travel time for a reflected wave between focus
        at depth Zf and station at depth Zs, separated by epicentral
        distance Repi.

        :param Zf:
        :param Zs:
        :param Repi:
        :param wave:
            see :meth:`find_reflection_ray_angles_and_tts`

        :return:
            (BL, tmin) tuple:
            - BL: int, bottom layer index
            - tmin: float, minimum travel time in seconds
        """
        if max(Zf, Zs) < self.max_depth:
            _, refl_tts = self.find_reflection_ray_angles_and_tts(
                Zf, Zs, Repi, wave=wave)
            BL = np.nanargmin(refl_tts)
            tmin = refl_tts[BL]
        else:
            BL, tmin = None, np.nan
        return (BL, tmin)

    def find_emerging_ray_angle_and_tt(self, Zf, Zs, Repi, wave='P'):
        """
        Determine angle and travel time of emerging ray (direct wave!) between
        focus at depth Zf and station at depth Zs, separated by epicentral
        distance Repi.

        :param Zf:
            float, depth of focus (in km)
        :param Zs:
            float, depth of station (in km)
        :param Repi:
            float, epicentral distance (= horizontal distance between
            earthquake focus and receiver station
        :param wave:
            see :meth:`get_velocities`

        Note: there is no constraint on the relative depth of Zf and Zs
        (they may be interchanged)

        :return:
            (theta, tt) tuple:
            - theta: float, angle of upgoing ray in radians
            - tt, float, travel time in seconds
        """
        # Determine higher/lower depth and layer index
        Zh = min(Zf, Zs)
        Zl = max(Zf, Zs)
        Lh = self.get_layer_index(Zh)
        Ll = self.get_layer_index(Zl)

        tt = None

        # TODO: if Zl coincides with layer boundary, consider
        # refracted wave traveling along this boundary
        # (not necessary, I think)

        if Lh == Ll:
            # Focus and station in same layer
            DZ = Zl - Zh
            theta = np.pi / 2 - np.arctan2(DZ, Repi)
            D = np.sqrt(DZ**2 + Repi**2)
            V = self.get_velocities(wave)
            tt = D / V[Lh]
        else:
            # Construct secondary velocity model going from Zh to Zl
            vmodel2 = self.constrain_depth_range(Zh, Zl)
            # print(vmodel2.depths)

            def minimize_func(theta):
                layer_angles = vmodel2.calc_layer_angles_above(
                    vmodel2.num_layers-1, theta, wave=wave)
                hdistances = vmodel2.calc_horizontal_distances(layer_angles)
                tot_hdistance = np.sum(hdistances)
                # print(tot_hdistance)
                return np.abs(Repi - tot_hdistance)

            result = minimize_scalar(
                minimize_func, bounds=(0, np.pi/2.), method='bounded')
            if not result.success:
                raise RuntimeError(result.message)
            theta = result.x
            layer_angles = vmodel2.calc_layer_angles_above(
                vmodel2.num_layers-1, theta, wave=wave)
            travel_times = vmodel2.calc_travel_times(
                layer_angles, wave=wave)
            tt = np.sum(travel_times)
        if tt is None:
            raise RuntimeError(
                'Could not determine emerging ray angle and travel time '
                'for direct wave.'
            )
        return (theta, tt)

    def calc_dirwav_tt(self, Zf, Zs, Repi, wave='P'):
        """
        Compute travel time of direct wave between focus at depth Zf
        and station at depth Zs, separated by epicentral distance Repi.

        :param Zf:
        :param Zs:
        :param Repi:
        :param wave:
            see :meth:`find_emerging_ray_angle_and_tt`

        Note: there is no constraint on the relative depth of Zf and Zs
        (they may be interchanged)

        :return:
            float, travel time in seconds
        """
        return self.find_emerging_ray_angle_and_tt(Zf, Zs, Repi, wave=wave)[1]

    def calc_all_tt(self, Zf, Zs, Repi):
        """
        Compute travel times for all simple phases:
        - Pg: direct P wave
        - Pn: refracted P wave
        - Pm: reflected P wave
        - P: min(Pg, Pn, Pm)
        - Sg: direct S wave
        - Sn: refracted S wave
        - Sm: reflected S wave
        - S: min(Sg, Sn, Sm)

        :param Zf:
        :param Zs:
        :param Repi:
            see :meth:`calc_dirwav_tt` or :meth:`calc_min_refwav_tt`

        :return:
            dict, mapping phase names (str) to travel times in seconds (float)
        """
        phase_tt = {
            'Pg': self.calc_dirwav_tt(Zf, Zs, Repi, wave='P'),
            'Pn': self.calc_min_refwav_tt(Zf, Zs, Repi, wave='P')[1],
            'Pm': self.calc_min_reflection_tt(Zf, Zs, Repi, wave='P')[1],
            'Sg': self.calc_dirwav_tt(Zf, Zs, Repi, wave='S'),
            'Sn': self.calc_min_refwav_tt(Zf, Zs, Repi, wave='S')[1],
            'Sm': self.calc_min_reflection_tt(Zf, Zs, Repi, wave='S')[1],
        }
        phase_tt['P'] = min(phase_tt['Pg'], phase_tt['Pn'], phase_tt['Pm'])
        phase_tt['S'] = min(phase_tt['Sg'], phase_tt['Sn'], phase_tt['Sm'])
        return phase_tt

    def calc_min_tt(self, Zf, Zs, Repi, wave='P'):
        """
        Compute minimum travel time (refracted or direct wave) between
        focus at depth Zf and station at depth Zs, separated by
        epicentral distance Repi.

        :param Zf:
        :param Zs:
        :param Repi:
        :param wave:
            see :meth:`calc_dirwav_tt` or :meth:`calc_min_refwav_tt`

        :return:
            (tmin, wave_type) tuple:
            - tmin: float, minimum travel time in seconds
            - wave_type: str, type of wave ('REF' or 'DIR')
        """
        ref_tmin = self.calc_min_refwav_tt(Zf, Zs, Repi, wave=wave)[1]
        dir_tmin = self.calc_dirwav_tt(Zf, Zs, Repi, wave=wave)
        return (ref_tmin, 'REF') if ref_tmin < dir_tmin else (dir_tmin, 'DIR')

    def calc_takeoff_and_incidence_angles(self, Zf, Zs, Repi, wave='P',
                                          wave_type=None):
        """
        Compute takeoff and incidence angles between
        focus at depth Zf and station at depth Zs, separated by
        epicentral distance Repi.

        :param Zf:
        :param Zs:
        :param Repi:
        :param wave:
            see :meth:`calc_dirwav_tt` or :meth:`calc_min_refwav_tt`
        :param wave_type:
            str, type of wave ('REF', 'DIR', 'REFL' or 'N', 'G', 'M')
            Alternatively, wave type may also be included in :param:`wave`,
            e.g. as 'Pg', 'Pn', 'Pm'
            (default: None, will calculate for fastest arrival)

        :return:
            (takeoff_angle, incidence_angle) tuple of floats:
            - takeoff_angle: angle at which ray leaves the source,
                measured from downwards pointing vertical (in degrees)
            - incidence_angle, angle at which ray arrives at the station,
                measured from downards pointing vertical (in degrees)
        """
        # Not the most efficient way, some calculations are probably duplicated
        if len(wave) > 1 and wave_type is None:
            wave, wave_type = wave[:1], wave[1:].upper()

        if not wave_type:
            wave_type = self.calc_min_tt(Zf, Zs, Repi, wave=wave)[1]
        wave_type = wave_type.upper()

        if wave_type in ('REF', 'N'):
            BL, _tmin = self.calc_min_refwav_tt(Zf, Zs, Repi, wave=wave)
            if BL is None:
                return (None, None)
            takeoff_angle = self.calc_critical_angles(wave)[BL]
        elif wave_type in ('DIR', 'G'):
            BL = self.get_layer_index(max(Zf, Zs))
            takeoff_angle = self.find_emerging_ray_angle_and_tt(
                Zf, Zs, Repi, wave=wave)[0]
        elif wave_type in ('REFL', 'M'):
            BL, _tmin = self.calc_min_reflection_tt(Zf, Zs, Repi, wave=wave)
            takeoff_angle = self.find_reflection_ray_angles_and_tts(
                Zf, Zs, Repi, wave=wave)[0][BL]
        else:
            raise ValueError(f'Unknown wave type: {wave_type}')
        incidence_angle = self.calc_layer_angles_above(
            BL, takeoff_angle, wave=wave)[0]
        # Next line converts to angle from horizontal downwards,
        # but this is not the convention
        # (see https://service.iris.edu/irisws/rotation/docs/1/help/)
        # incidence_angle = np.pi/2. - incidence_angle

        # Direct wave leaves upwards from source
        if wave_type in ('DIR', 'G'):
            takeoff_angle = np.pi - takeoff_angle

        return np.degrees(takeoff_angle), np.degrees(incidence_angle)

    def calc_tt_and_angles(self, Zf, Zs, Repi, wave='P', wave_type=None):
        """
        Simultaneously compute travel time and takeoff and incidence angles
        between focus at depth Zf and station at depth Zs, separated by
        epicentral distance Repi.

        This is more efficient than calling :math:`calc_min_tt` and
        :meth:`calc_takeoff_and_incidence_angles` separately

        :param Zf:
        :param Zs:
        :param Repi:
        :param wave:
            see :meth:`calc_dirwav_tt` or :meth:`calc_min_refwav_tt`
        :param wave_type:
            str, type of wave ('REF', 'DIR', 'REFL' or 'N', 'G', 'M')
            Alternatively, wave type may also be included in :param:`wave`,
            e.g. as 'Pg', 'Pn', 'Pm'
            (default: None, will calculate for fastest arrival)

        :return:
            (travel_time takeoff_angle, incidence_angle, wave_type) tuple:
            - travel_time: float, travel time (in seconds)
            - takeoff_angle: float, angle at which ray leaves the source,
                measured from downwards pointing vertical (in degrees)
            - incidence_angle, float, angle at which ray arrives at
                the station, measured from downwards pointing vertical
                (in degrees)
            - wave_type: str, corresponding wave type
        """
        if len(wave) > 1 and wave_type is None:
            wave, wave_type = wave[:1], wave[1:].upper()

        # Use specified wave type or find fastest
        wave_types = [wave_type] if wave_type else ['REF', 'DIR']

        travel_times = []
        takeoff_angles = []
        incidence_angles = []
        for wave_type in wave_types:
            if wave_type in ('DIR', 'G'):
                to_angle, tt = self.find_emerging_ray_angle_and_tt(
                    Zf, Zs, Repi, wave=wave)
                # Direct wave leaves upwards from source
                to_angle = np.pi - to_angle
                BL = self.get_layer_index(max(Zf, Zs))
            elif wave_type in ('REF', 'N'):
                BL, tt = self.calc_min_refwav_tt(Zf, Zs, Repi, wave=wave)
                if BL is None:
                    to_angle, tt = np.nan, np.nan
                else:
                    to_angle = self.calc_critical_angles(wave)[BL]
            elif wave_type in ('REFL', 'M'):
                if max(Zf, Zs) < self.max_depth:
                    _to_angles, _tts =\
                        self.find_reflection_ray_angles_and_tts(
                            Zf, Zs, Repi, wave=wave)
                    BL = np.nanargmin(_tts)
                    tt = _tts[BL]
                    to_angle = _to_angles[BL]
                else:
                    tt, to_angle = np.nan, np.nan
                    BL = None
            else:
                raise ValueError(f'Unknown wave type: {wave_type}')
            travel_times.append(tt)
            takeoff_angles.append(to_angle)

            # Calculate incidence angle from takeoff angle and bottom layer idx
            if BL is not None:
                incidence_angle = self.calc_layer_angles_above(
                    BL, to_angle, wave=wave)[0]
            else:
                incidence_angle = np.nan
            incidence_angles.append(incidence_angle)

        takeoff_angles = np.degrees(takeoff_angles)
        incidence_angles = np.degrees(incidence_angles)

        # Find index corresponding to fastest arrival
        # or use the only one available
        i = np.nanargmin(travel_times) if len(travel_times) > 1 else 0

        return (
            travel_times[i], takeoff_angles[i], incidence_angles[i],
            wave_types[i]
        )

    def calc_sp_interval(self, Zf, Zs, Repi):
        """
        Compute time difference between S and P-wave arrivals

        :param Zf:
        :param Zs:
        :param Repi:
            see :meth:`calc_dirwav_tt` or :meth:`calc_min_refwav_tt`

        :return:
            float, S-P time interval
        """
        Pmin_tt = self.calc_min_tt(Zf, Zs, Repi, wave='P')[0]
        Smin_tt = self.calc_min_tt(Zf, Zs, Repi, wave='S')[0]
        return Smin_tt - Pmin_tt

    def calc_epicentral_distance(self, Zf, Zs, sp_interval):
        """
        Compute epicentral distance based on difference between
        S-wave and P-wave arrival times

        :param Zf:
            float, depth of focus (in km)
        :param Zs:
            float, depth of station (in km)
        :param sp_interval:
            float, time between S-wave and P-wave arrivals (in seconds)

        :return:
            float, epicentral distance (in km)
        """
        def minimize_func(repi):
            SminusP = self.calc_sp_interval(Zf, Zs, repi)
            return np.abs(SminusP - sp_interval)

        max_dist = sp_interval * 10
        result = minimize_scalar(
            minimize_func, bounds=(0, max_dist), method='bounded')
        if result.success:
            Repi = result.x
        else:
            raise RuntimeError(result.message)
        return Repi

    def calc_tt_residuals(self, Zf, Tf, Repi, Zs, Ts, Vidx):
        """
        Compute travel-time residuals.
        Residuals are positive when observed travel time is larger
        than calculated travel time

        :param Zf:
            float, depth of focus (in km)
        :param Tf:
            float, origin time (in seconds)
        :param Repi:
            1D array, epicentral distances (in km) [num_phases]
        :param Zs:
            1D array, station depths (in km) [num_phases]
        :param Ts:
            1D array, phase arrival times (in seconds) [num_phases]
        :param Vidx:
            1D array, velocity indexes (0=P or 1=S) [num_phases]

        :return:
            (tt, tt_residuals) tuple
            - tt: 1D array, calculated travel times (in seconds)
            - tt_residuals: 1D array, travel time residuals (in seconds)
        """
        assert len(Repi) == len(Zs) == len(Ts) == len(Vidx)

        num_phases = len(Repi)
        tt = np.zeros(num_phases)
        tt_residuals = np.zeros(num_phases)
        for k in range(num_phases):
            repi, zs, ts, vidx = Repi[k], Zs[k], Ts[k], Vidx[k]
            wave = 'PS'[vidx]
            tt_obs = ts - Tf
            tt_calc = self.calc_min_tt(Zf, zs, repi, wave=wave)[0]
            tt[k] = tt_calc
            tt_residuals[k] = tt_obs - tt_calc

        return tt, tt_residuals

    def calc_path_elements(self, Zf, Zs, Repi, wave_type, wave='P'):
        """
        Compute elements (X and Y coordinates) of travel path

        :param Zf:
        :param Zs:
        :param Repi:
            see :meth:`calc_min_tt`
        :param wave_type:
            str, wave type: 'DIR', 'REF' or 'REFL' for direct, refracted
            and reflected wave
        :param wave:
            see :meth:`calc_min_tt`

        :return:
            X, Y tuple of arrays: X and Y coordinates (in km) of travel path,
            starting from the hypocenter (0, Zf)
        """
        X, Y = [], []

        if wave_type in ('DIR', 'g'):
            # Direct wave
            X, Y = [Repi], [Zs]
            if not np.isclose(Zs, Zf):
                vmodel2 = self.constrain_depth_range(Zs, Zf)
                bottom_angle, _tt = vmodel2.find_emerging_ray_angle_and_tt(
                    vmodel2.max_depth, 0, Repi, wave=wave)
                layer_angles = vmodel2.calc_layer_angles_above(
                    vmodel2.num_layers - 1, bottom_angle, wave=wave)
                hdistances = vmodel2.calc_horizontal_distances(layer_angles)
                for layer, hdist in enumerate(hdistances):
                    X.append(X[-1] - hdist)
                    Y.append(vmodel2.depths[layer+1] + Zs)
            else:
                X.append(0)
                Y.append(Zf)
            X = X[::-1]
            Y = Y[::-1]

        elif wave_type in ('REF', 'REFR', 'n'):
            # Refracted wave
            BL, _tmin = self.calc_min_refwav_tt(Zf, Zs, Repi, wave=wave)
            if BL is not None:
                bottom_angle = self.calc_critical_angles(wave=wave)[BL]
                vmodel1 = self.constrain_depth_range(Zf, self.depths[BL+1])
                layer_angles = vmodel1.calc_layer_angles_above(
                    vmodel1.num_layers - 1, bottom_angle, wave=wave)
                hdistances = vmodel1.calc_horizontal_distances(layer_angles)
                X.append(0)
                Y.append(Zf)
                for layer, hdist in enumerate(hdistances):
                    X.append(X[-1] + hdist)
                    Y.append(vmodel1.depths[layer+1] + Zf)
                vmodel2 = self.constrain_depth_range(Zs, self.depths[BL+1])
                layer_angles = vmodel2.calc_layer_angles_above(
                    vmodel2.num_layers - 1, bottom_angle, wave=wave)
                hdistances = vmodel2.calc_horizontal_distances(layer_angles)
                X.append(Repi - np.sum(hdistances))
                Y.append(Y[-1])
                for layer, hdist in enumerate(hdistances[::-1]):
                    X.append(X[-1] + hdist)
                    Y.append(
                        vmodel2.depths[vmodel2.num_layers - layer - 1] + Zs)

        elif wave_type in ('REFL', 'm'):
            # Reflected wave
            BL, _tmin = self.calc_min_reflection_tt(Zf, Zs, Repi, wave=wave)
            if BL is not None:
                refl_angles, _ = self.find_reflection_ray_angles_and_tts(
                    Zf, Zs, Repi, wave=wave)
                bottom_angle = refl_angles[BL]
                vmodel1 = self.constrain_depth_range(Zf, self.depths[BL+1])
                layer_angles = vmodel1.calc_layer_angles_above(
                    vmodel1.num_layers - 1, bottom_angle, wave=wave)
                hdistances = vmodel1.calc_horizontal_distances(layer_angles)
                X.append(0)
                Y.append(Zf)
                for layer, hdist in enumerate(hdistances):
                    X.append(X[-1] + hdist)
                    Y.append(vmodel1.depths[layer+1] + Zf)
                vmodel2 = self.constrain_depth_range(Zs, self.depths[BL+1])
                layer_angles = vmodel2.calc_layer_angles_above(
                    vmodel2.num_layers - 1, bottom_angle, wave=wave)
                hdistances = vmodel2.calc_horizontal_distances(layer_angles)
                for layer, hdist in enumerate(hdistances[::-1]):
                    X.append(X[-1] + hdist)
                    Y.append(
                        vmodel2.depths[vmodel2.num_layers - layer - 1] + Zs)

        else:
            raise ValueError(f'Unknown wave type: {wave_type}')

        return (np.array(X), np.array(Y))

    def calc_travel_distance(self, Zf, Zs, Repi, wave_type, wave='P'):
        """
        Compute travel distance

        :param Zf:
        :param Zs:
        :param Repi:
        :param wave_type:
        :param wave:
            see :meth:`calc_path_elements`

        :return:
            float, travel distance (in km)
        """
        X, Y = self.calc_path_elements(Zf, Zs, Repi, wave_type, wave=wave)
        return np.sum(np.hypot(np.diff(X), np.diff(Y)))

    def calc_average_velocity(self, Zf, Zs, Repi, wave='P'):
        """
        Compute average velocity over travel path

        :param Zf:
        :param Zs:
        :param Repi:
        :param wave:
            see :meth:`calc_min_tt`

        :return:
            float, average velocity (in km/s)
        """
        tt, wave_type = self.calc_min_tt(Zf, Zs, Repi, wave=wave)
        d = self.calc_travel_distance(Zf, Zs, Repi, wave_type, wave=wave)
        return d / tt
