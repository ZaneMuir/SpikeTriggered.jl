# SpikeTriggered.jl

Package for spike triggered neural characterization.

## structure

- PSTH:
  - raster
  - psth2timestamps
  - TODO: psth2trace
  - histogram_gsl
  - XXX: gsl_gaussian_filter
  - XXX: spk_filter
  - TODO: spike_gaussian_filter
- Stats:
  - TODO: spike_triggered_average
  - TODO: spike_triggered_covariance
  - spike_train_spectrum_power
- OnlineSpikes
  - detect_cross_pnt
  - detect_events
- TODO: Waveforms
  - TODO: trace_event_related
  - TODO: trace_downsample
  - TODO: trace_upsample
  - TODO: trace_gaussian_smooth
  - TODO: trace_remove_spike
  - TODO: trace_remove_event
  - TODO: make_event_template
- Stimulus
  - SparseNoise
