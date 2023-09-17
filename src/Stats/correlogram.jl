
# @doc raw"""
#     correlogram(target::Vector{T}, reference::Vector{T}; bin_size, window_size, is_auto) where {T} -> (Vector, Vector)

# calculate the correlogram between the target and reference timestamps.

# ## Arguments:
# - `target::Vector{T}`, `reference::Vector{T}`: timestamps of events (NB: the calculation is based on the timestamps, not the histograms.)

# ## Keyword Arguments:
# - `bin_size`: bin size in second [default: 1e-2]
# - `window_size`: total number of bins [default: 301]
# - `is_auto`: flag for autocorrelogram, i.e. ignoring the same timestamp. [default: false]

# One should use either two quick access functions: `autocorrelogram` and `crosscorrelogram`.

# ## References
# - [wikipedia](https://en.wikipedia.org/wiki/Correlogram)
# - [analysis tools for MULAB wabpage](https://www.med.upenn.edu/mulab/analysis.html#Introduction)
# """
# function correlogram(target::Vector{T}, reference::Vector{T}; bin_size=1e-2, window_size=301, is_auto=false) where {T}
#     _K_half = div(window_size, 2)
#     _bin_boundary = (-_K_half:_K_half) .* bin_size  #FIXME: boundary not right
#     output = zeros(Int64, window_size-1)
#     for pt in reference
#         _offset = target .- pt
#         _offset = is_auto ? _offset[_offset .!= 0] : _offset
#         output .+= histogram_gsl(Float64.(_offset), Float64.(_bin_boundary))
#     end
#     _bin_boundary[2:end] .- bin_size/2, output
# end

# autocorrelogram(ts_vec; kwargs...) = correlogram(ts_vec, ts_vec; is_auto=true, kwargs...)
# crosscorrelogram(args...; kwargs...) = correlogram(args...; is_auto=false, kwargs...)

function cross_correlation(f::AbstractArray{T}, g::AbstractArray{T}, dims=1) where {T <: Real}
    N = size(f, dims)
    @assert reduce(&, size(f) .== size(g)) "input sizes are not matched."
    fft_plan = FFTW.plan_rfft(f, dims)

    F_bar = fft_plan * f .|> conj
    G = fft_plan * g

    _xcorr = abs.(FFTW.irfft(F_bar .* G, N, dims))
    FFTW.fftshift(_xcorr, dims)
end

function cross_correlation_window(f::AbstractArray, dims=1; fs=1)
    N = size(f, dims)
    _t = FFTW.fftfreq(N, N / fs)
    FFTW.fftshift(_t)
end

function auto_correlation(f, dims=1; corrected=false)
    _corr = cross_correlation(f, f, dims)
    if corrected
        _corr = FFTW.ifftshift(_corr)
        _corr[1] -= sum(_corr)
        FFTW.fftshift(_corr)
    else
        _corr
    end
end


@doc raw"""

The covariogram is defined as:

```math
V \equiv \langle (S_1^{r} - P_1) \odot (S_2^{r} - P_2) \rangle = \langle S_1^{r} \odot S_2^{r} \rangle - P_1 \odot P_2
```
"""
function covariogram(S1::AbstractMatrix{T}, S2::AbstractMatrix{T}; fs=1) where {T <: Real}
    N = size(S1, 2)
    P1 = mean(S1, dims=2)[:]
    P2 = mean(S2, dims=2)[:]
    σsq1 = var(S1, dims=2)[:]
    σsq2 = var(S2, dims=2)[:]

    raw_xcorr_mat = cross_correlation(S1, S2, 1)
    raw_xcorr = mean(raw_xcorr_mat, dims=2)[:]
    shuffle_corr = cross_correlation(P1, P2)
    V = raw_xcorr .- shuffle_corr
    σsqV = (cross_correlation(σsq1, σsq2) .+ cross_correlation(P1.^2, σsq2) .+ cross_correlation(σsq1, P2.^2)) ./ N
    σV = sqrt.(σsqV)
    trange = cross_correlation_window(S1, 1; fs)
    (; t=trange, V, σ=σV)
end