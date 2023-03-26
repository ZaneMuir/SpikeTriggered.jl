@doc raw"""
    spike_filter(spk::Vector{T}, proj, kernel::Function; norm::Bool=true) where {T <: Real} -> Vector{T}

generating smoothed curve from spike trains. equivalent to the smoothed PSTH.

## Arguments:
- spk: spike train vector
- proj: range or vector of timestamps interested
- kernel: function of the smoothing kernel, in format of `T -> T`
"""
function spike_filter(spk::Vector{T}, proj, kernel::Function; norm::Bool=true) where {T <: Real}
    isempty(spk) && (return zeros(T, size(proj)))
    
    _psth = zeros(T, length(proj))
    for idx in eachindex(proj)
        @inbounds @fastmath _psth[idx] = sum(kernel, spk .- proj[idx])
    end
    if norm
        _psth ./ length(spk)
    else
        _psth
    end
end

@deprecate spk_filter(spk, proj, kernel; norm) spike_filter(spk, proj, kernel; norm)

@doc raw"""
    spike_gaussian_filter(spk::Vector{T}, proj; σ=0.03) where {T<:Real} -> Vector{T}

Gaussian smoothing for a given spike train.

```math
\text{G}(; \sigma) = \frac{1}{\sqrt{2 \pi} \sigma} \exp(- \frac{x^2}{2 \sigma^2})
```

## Arguments
- `spk::Vector{T}`: spike train vector.
- `proj`: range or vector of timestamps.

## Keyword Arguments:
- `σ`: sigma of the gaussian kernel (in seconds) [default: 0.03]

## Returns:
- `Vector{T}`: smoothed psth with the same length as `proj`.
"""
function spike_gaussian_filter(spk::Vector{T}, proj; σ=0.005, kwargs...) where {T <: Real}
    k_gaussian = (x::T) -> exp(x^2 / (- 2 * σ ^ 2)) / (σ * sqrt(2 * π))
    
    spike_filter(spk, proj, k_gaussian; norm=false, kwargs...)
end