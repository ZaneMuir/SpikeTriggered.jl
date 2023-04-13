# SpikeTriggered.jl Documentation

The [SpikeTriggered](https://github.com/ZaneMuir/SpikeTriggered.jl) package provides a collection of functions for analysis related to spike trains. Including:

- converting spike train into rasters.
- getting PSTH from spike train, by using histogram or Gaussian filter.
- estimating RFs using forward correlation and reserve correlation.
- thresholding waveforms to extract spike events.
- etc.

Functions for spike train statistics are in the `Stats` submodule and functions for waveform processing are in the `Waveforms` submodule. Functions in the `Stimulus` submodule are for internal uses of HirschLab.

For visualization, checkout the [`NeuroMakie.jl`](https://github.com/ZaneMuir/NeuroMakie.jl) package.
