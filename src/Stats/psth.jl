
#### Histogram

@doc raw"""
    histogram_gsl(u_arr, edges) where {T <: Real}

Count histogram by providing the edges. If provides `n+1` edges,
it would return `n` length histogram.
Edges would work as: left bound <= value < right bound.

# Arguments
- `u_arr`: a vector of values
- `edges`: a vector of edges

Note: StatsBase.Histogram is slower than GSL but with less memory footprint.
"""
function histogram_gsl(u_arr::AbstractVector{T}, edges) where {T <: Real}
    u_arr = Cdouble.(u_arr)
    edges = Cdouble.(edges)

    edges = sort(edges)

    n = length(edges)-1
    gsl_hist = GSL.histogram_alloc(n)
    GSL.histogram_set_ranges(gsl_hist, edges, n+1)

    for idx in eachindex(u_arr)
        @inbounds GSL.histogram_increment(gsl_hist, u_arr[idx])
    end

    myhist = zeros(Int, n)
    for idx in 1:n
        @inbounds myhist[idx] = GSL.histogram_get(gsl_hist, idx-1)
    end

    GSL.histogram_free(gsl_hist)

    return myhist
end

@doc raw"""
    spike_histogram(spk::AbstractVector, edges::AbstractVector)

PSTH as histogram, where edges are defined by user.

# Arguments:

- `spk::AbstractVector{T}`: spike train
- `edges::AbstractVector`: edges of histogram, noted that it should be length of `n+1`.

# Returns:

- psth as `Vector{T}` of length `n`.
"""
function spike_histogram(spk::AbstractVector{T}, edges::AbstractVector) where {T <: Real}
    if isempty(spk)
        zeros(Int, length(edges)-1)
    else
        histogram_gsl(spk, edges)
    end
    #not so fast histogram: _hist = StatsBase.fit(StatsBase.Histogram, spk, edges, closed=:left).weights
end

@doc raw"""
    spike_histogram(spk::AbstractVector, binsize::Real)

PSTH as histogram, where edges are created based on the provided binsize.

# Arguments:

- `spk::AbstractVector{T}`: spike train
- `binsize::Real`: binsize of histogram. edges will be created as: `range(_min_time, step=binsize, length=_L+1)`
"""
function spike_histogram(spk::AbstractVector{T}, binsize::Real) where {T <: Real}
    (_min, _max) = extrema(spk)
    _L = floor(Int, (_max - _min) / binsize)
    _edges = range(_min, step=binsize, length=_L+1)
    spike_histogram(spk, _edges)
end

@doc raw"""
    spike_histogram(raster::Vector{Vector{T}}, args...; norm=true)

PSTH as histogram from averaging trials.
"""
function spike_histogram(raster::Vector{Vector{T}}, args...; norm::Bool=true) where {T <: Real}
    _flatten = reduce(vcat, raster; init=T[])
    _psth = spike_histogram(_flatten, args...)
    norm ? _psth ./ length(raster) : _psth
end

#### Filter

@doc raw"""
    spike_filter(spk, proj, kernel::Function; kwargs...) where {T <: Real} -> Vector{T}

generating smoothed curve from spike trains. equivalent to convolution.

```math
\text{PSTH}(t) = (h * s)(t) = \sum_i \delta(t - t'_i) h (t - t'_i),\; t \in \text{proj}
```

# Arguments:

- spk: spike train vector
- proj: range or vector of timestamps interested
- kernel: function of the smoothing kernel, should be `(::T; kwargs...)::T`

# Returns:

- psth as `Vector{T}`, same length as `proj`.
"""
function spike_filter(spk::AbstractVector{T}, proj::AbstractVector, kernel::Function; kwargs...) where {T <: Real}

    isempty(spk) && (return zeros(T, size(proj)))

    # sacrifice memory for speed
    _psth = sum(kernel.(T.(proj)' .- spk; kwargs...); dims=1)[:]
end

@doc raw"""
    spike_filter(raster, proj, kernel::Function; kwargs...) where {T <: Real} -> Vector{T}

generating smoothed curve from rasters.
"""
function spike_filter(raster::Vector{Vector{T}}, args...; kwargs...) where {T <: Real}
    _flatten = reduce(vcat, raster)
    _N = length(raster)
    _psth = spike_filter(_flatten, args...; kwargs...) ./ _N
end

@doc raw"""
    gaussian_kernel(x::T; σ::T)

Inline [Gaussian function](https://en.wikipedia.org/wiki/Gaussian_filter).

```math
g(x) = \frac{\exp(-\frac{x^2}{2 σ^2})}{\sqrt{2\pi} σ}
```
"""
@inline gaussian_kernel(x::T; σ::T) where{T <: AbstractFloat} = exp(x^2 / (- 2 * σ ^ 2)) / (σ * sqrt(2 * T(π)))

@doc raw"""
    spike_filter_gaussian(spk_or_raster, proj; σ=0.010)

Spike filter using [Gaussian function](https://en.wikipedia.org/wiki/Gaussian_filter).

```math
g(x) = \frac{\exp(-\frac{x^2}{2 σ^2})}{\sqrt{2\pi} σ}
```
"""
function spike_filter_gaussian(spk_or_raster::Union{AbstractVector{T}, Vector{Vector{T}}}, proj::AbstractVector; σ::Real=0.010) where {T <: Real}
    spike_filter(spk_or_raster, proj, gaussian_kernel; σ=T(σ))
end

@doc raw"""
    get_histogram_center(edges)

get the center of edge vector. make sure edges are sorted.
"""
function get_histogram_center(edges::AbstractVector)
    edges[1:end-1] .+ diff(edges) ./ 2
end

@doc raw"""
    spike_xcorr(target, reference, roi) where {T <: Real}

The crosscorrelogram shows a count of the spikes of the target cell
at specific time delays with respect the spikes of the reference cell.

If there are multiple trials, use the `spike_xcorr_shifted` to correct
any bias with a shift predictor.
"""
function spike_xcorr(target::AbstractVector{T}, reference::AbstractVector{T}, roi::AbstractVector) where {T <: Real}
    rez = zeros(Int, length(roi)-1)
    for ref_t in reference
        _my_psth = spike_histogram(target .- ref_t, roi)
        rez .+= _my_psth
    end
    rez
end

@doc raw"""
    spike_xcorr_shifted(target, reference, roi; shift_t, shift_t_end) where {T <: Real}

The crosscorrelogram shows a count of the spikes of the target cell
at specific time delays with respect the spikes of the reference cell.
Corrected with a shift predictor.

## Keyword Arguments:
- shift_t: Δt of each shift, usually shift the spike by trials
- shift_t_end: the end limit of spike_time, usually the end time of the last trial. [default: maximum(reference)]
"""
function spike_xcorr_shifted(target::AbstractVector{T}, reference::AbstractVector{T}, roi::AbstractVector;
        shift_t::AbstractVector{T},
        shift_t_end::Union{Nothing, T}=nothing
        ) where {T <: Real}
    raw = spike_xcorr(target, reference, roi)
    shift_t_end = isnothing(shift_t_end) ? maximum(reference) : shift_t_end
    shift_predictor_mat = zeros(Int, length(roi)-1, length(shift_t))
    for (idx, shift_it) in enumerate(shift_t)
        _shifted_ref = reference .- shift_it
        _shifted_ref[_shifted_ref .< 0] .+= .+ shift_t_end
        shift_predictor_mat[:, idx] = spike_xcorr(target, _shifted_ref, roi)
    end
    raw .- mean(shift_predictor_mat, dims=2)[:], shift_predictor_mat
end

#TODO: crosscorrelogram from rasters, spike_xcorr_shifted(target_raster, reference_raster, roi::AbstractVector)
