# Raster and PSTH

Raster plot and [PSTH](https://en.wikipedia.org/wiki/Peristimulus_time_histogram) are very useful for neural activities to periodic stimuli. Patterns of firing can be easily seen and further analyized when spikes are aligned with the onset of the stimuli.

## Getting Rasters

Rasters will be created as `Vector{Vector{T}}`, where each subvector will be spike times relative to the stimulus onset.
Besides the spike train, you should at least have the onset time of the stimulus (`marker::AbstractVector`), and specify the length of time before the onset and duration of the stimulus and time after the stimulus.

```julia
julia> spike_train = rand(1000) * 50
julia> marker = 0:10:40
julia> raster = SpikeTriggered.Stats.spike_raster(spike_train, marker;
                                           head=0.5, duration=1.0, tail=0.5)
...
```

## Getting PSTH

### Using Histogram

Just like get a histogram from any other data. You could either define the bin size of your histogram:

```julia
julia> bin_size = 0.5
julia> psth = SpikeTriggered.Stats.spike_histogram(spike_train/raster, bin_size)
```

Or providing the edges of your bins:

```julia
julia> edges = 0:0.01:50.0
julia> psth = SpikeTriggered.Stats.spike_histogram(spike_train/raster, edges)
```

The latter is preferred, as you could have finer control and you could have a nonlinear timescale. Just don't forget to append the last right edge.

### Using Filter

Sometimes, it's arbitary to find a good bin size, or you want to smooth the histogram, you can use a filter to generate a smooth PSTH. Essentially, it is a convolution over time steps provided by `projection::AbstractVector`.

```math
\text{PSTH}(t) = (h * s)(t) = \sum_i \delta(t - t'_i) h (t - t'_i),\; t \in \text{projection}
```

where `h` is the kernel function, `t'_i` is the spike time.

Gaussian kernel is mostly common one, so you can generate it by:

```julia
julia> psth = SpikeTriggered.Stats.spike_filter_gaussian(spike_train/raster, projection; σ=0.010)
```

Or, if you want, you can also create your own kernel:

```julia
julia> my_kernel(dt; thresh) = dt < thresh ? 0 : 1
julia> psth = SpikeTriggered.Stats.spike_filter(spike_train/raster, projection, my_kernel; thresh=0.5)
```

## Functions

```@docs
SpikeTriggered.Stats.spike_raster
SpikeTriggered.Stats.spike_histogram
SpikeTriggered.Stats.spike_filter
SpikeTriggered.Stats.spike_filter_gaussian
SpikeTriggered.PSTH.psth2timestamps
SpikeTriggered.PSTH.psth2trace
```
