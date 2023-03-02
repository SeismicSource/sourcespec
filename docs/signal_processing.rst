.. _signal_processing:

#################
Signal Processing
#################

The following documentation explains, in chronological order, all the steps
performed to construct the amplitude spectra used for the inversion.

Trace Processing
~~~~~~~~~~~~~~~~

1. Traces are checked for gaps and overlaps. Traces with cumulative gap or
   overlap duration larger than ``gap_max``  or ``overlap_max``, respectively,
   are skipped.

2. The trace mean is removed (only if the configuration parameter
   ``trace_units`` is set to ``auto``, i.e., if the trace has not been
   preprocessed outside SourceSpec).

3. Gaps and overlaps are merged.

4. Traces with RMS smaller than ``rmsmin`` are skipped.

5. Traces are optionally checked for clipping (see :ref:`clipping_detection`).

6. Instrumental response is removed and trace transformed in its physical units
   (e.g., velocity, acceleration).

7. Trace is filtered.

8. Signal to noise ratio is measured as the ratio between signal RMS in the P-
   or S-window (depending on the ``wave_type`` parameter) and the RMS of the
   noise window. Traces with signal to noise ratio smaller than ``sn_min`` are
   skipped.

See the source code of :meth:`ssp_process_traces.process_traces` for
implementation details.

Spectral Processing
~~~~~~~~~~~~~~~~~~~

1. A check is made to see if the signal window contains enough data. Traces
   with no data or more than 25% of zero values are skipped. This is typically
   due to a wrong specification of signal window (e.g., wrong P or S arrivals).

2. Traces are integrated to displacement in time domain, if ``time_domain_int``
   is ``True``.

3. Signal and noise windows are cut and tapered.

4. The noise window is checked. If the ratio between noise RMS and signal RMS
   is smaller than 1e-6, the noise is considered not significant and the trace
   is skipped.

5. The amplitude spectra of the signal and noise windows are computed, using
   :func:`numpy.fft.rfft()`.
   If ``spectral_win_length`` is not ``None``, the signal is zero-padded to
   this length before computing the Fast Fourier Transform.

6. Spectra are integrated to displacement in frequency domain, if
   ``time_domain_int`` is ``False``.

7. Amplitude spectra are windowed (see config parameters ``freq1_broadb``,
   ``freq2_broadb`` and similar).

8. Geometrical spreading is corrected (see
   :ref:`theoretical_background:Geometrical Spreading`).

9. Amplitude spectra are converted to seismic moment units (see
   :ref:`theoretical_background:Building Spectra`).

10. Amplitude spectra are resampled in :math:`log_{10}` frequency spacing.

11. The resampled spectra are smoothed. The smoothing window width is defined
    in frequency decades (config parameter ``spectral_smooth_width_decades``)

12. The spectral signal to noise ratio is checked. Amplitude spectra with
    signal to noise ratio smaller than ``spectral_sn_min`` are skipped.

13. The "H" component is built based on one or more spectral components,
    depending on the ``wave_type`` and ``ignore_vertical`` config parameters:

    - if ``wave_type`` is ``S``:

        - if ``ignore_vertical`` is ``False``, the two horizontals and the
          vertical components are combined;
        - if ``ignore_vertical`` is ``True``, only the two horizontals
          components are combined;

    - if ``wave_type`` is ``SV``:

        - if ``ignore_vertical`` is ``False``, the radial and the vertical
          component are combined;
        - if ``ignore_vertical`` is ``True``, only the radial component is used;

    - if ``wave_type`` is ``SH``:

        - only the transverse component is used;

    - if ``wave_type`` is ``P``:

        - only the vertical component is used, independently from the value
          of ``ignore_vertical``.

    Spectra are combined through the root sum of squares (see
    :ref:`theoretical_background:Overview`).

14. All the amplitude spectra are converted to moment magnitude units (see
    :ref:`theoretical_background:Building Spectra`).

15. Station corrections are applied, if requested (see
    :ref:`theoretical_background:Station Residuals`).

16. The weight spectrum is built, depending on the config option ``weighting``.


See the source code of :meth:`ssp_build_spectra.build_spectra` for
implementation details.