# SpikeTriggered.jl

Package for spike triggered neural characterization.

## structure

- PSTH:
  - raster
  - psth2timestamps
  - psth2trace
  - histogram
  - histogram_gsl
  - spike_gaussian_filter
    - spk_filter
    - k_gaussian
- Stats:
  - spike_triggered_average
  - TODO: spike_triggered_covariance
  - spike_train_spectrum_power
  - detect_burst
  - interp_burst
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
