# SpikeTriggered.jl

Package for spike triggered neural characterization.

## Index:

- PSTH:
  - raster
  - psth2timestamps
  - psth2trace
  - histogram
  - histogram_gsl
  - spike_gaussian_filter
    - spike_filter
    - k_gaussian
- Stats:
  - spike_xcorr
  - spike_xcorr_shifted
  - spike_jpsth_raw
  - spike_jpsth_shifted
  - spike_triggered_average
  - spike_triggered_covariance
  - spike_train_spectrum_power (F1)
  - detect_burst
    - detect_burst_trn
    - detect_burst_lgn
    - split_tonic_burst
    - interp_burst
  - dispersion_index (Fano factor)
  - correlogram #FIXME
    - autocorrelogram
    - crosscorrelogram
- Waveforms
  - detect_cross_pnt
  - detect_events
  - gaussian_filter_gsl
  - filter_gaussian
  - filter_median
  - filter_LP_butterworth
  - filter_BP_butterworth
  - filter_HP_butterworth
  - trace_event_related
  - TODO: trace_downsample
  - TODO: trace_upsample
  - TODO: trace_remove_spike
  - TODO: trace_remove_event
  - TODO: make_event_template
- Stimulus
  - marker validation
    - process_first_of_paired_markers
    - remove_invalid_timestamps
    - remove_invalid_serial_timestamps!
    - remove_invalid_paired_timestamps!
  - apply_video_nonlinearity
  - SparseNoise
    - load_snra
    - parse_snra
    - Base.collect
  - TODO: GaussianNoise (for quick io loading)
