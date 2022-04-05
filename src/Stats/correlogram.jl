
@doc raw"""
    correlogram(target::Vector{T}, reference::Vector{T}; bin_size, window_size, is_auto) where {T} -> (Vector, Vector)

calculate the correlogram between the target and reference timestamps.

## Arguments:
- `target::Vector{T}`, `reference::Vector{T}`: timestamps of events (NB: the calculation is based on the timestamps, not the histograms.)

## Keyword Arguments:
- `bin_size`: bin size in second [default: 1e-2]
- `window_size`: total number of bins [default: 301]
- `is_auto`: flag for autocorrelogram, i.e. ignoring the same timestamp. [default: false]

One should use either two quick access functions: `autocorrelogram` and `crosscorrelogram`.

## References
- [wikipedia](https://en.wikipedia.org/wiki/Correlogram)
- [analysis tools for MULAB wabpage](https://www.med.upenn.edu/mulab/analysis.html#Introduction)
"""
function correlogram(target::Vector{T}, reference::Vector{T}; bin_size=1e-2, window_size=301, is_auto=false) where {T}
    _K_half = div(window_size, 2)
    _bin_boundary = (-_K_half:_K_half) .* bin_size  #FIXME: boundary not right
    output = zeros(Int64, window_size-1)
    for pt in reference
        _offset = target .- pt
        _offset = is_auto ? _offset[_offset .!= 0] : _offset
        output .+= histogram_gsl(Float64.(_offset), Float64.(_bin_boundary))
    end
    _bin_boundary[2:end] .- bin_size/2, output
end

autocorrelogram(ts_vec; kwargs...) = correlogram(ts_vec, ts_vec; is_auto=true, kwargs...)
crosscorrelogram(args...; kwargs...) = correlogram(args...; is_auto=false, kwargs...)