k_gaussian(σ::Real=1) = (x) -> exp(x^2 / (- 2 * σ ^ 2)) / (σ * sqrt(2 * π))

function spk_filter(spk::Vector{T}, proj, kernel::Function) where {T <: Real}
    _psth = zeros(T, length(proj))
    for idx in 1:length(proj)
        @inbounds @fastmath _psth[idx] = sum(kernel, spk .- proj[idx])
    end
    return _psth
end

@doc raw"""
    spike_gaussian_filter(spk::Vector{T}, proj; σ=0.03) where {T<:Real} -> Vector{T}

Gaussian smoothing for a given spike train.

## Arguments
- `spk::Vector{T}`: spike train vector.
- `proj`: range or vector of timestamps.

## Keyword Arguments:
- `σ`: sigma of the gaussian kernel (in seconds) [default: 0.03]

## Returns:
- `Vector{T}`: smoothed psth with the same length as `proj`.
"""
function spike_gaussian_filter(spk::Vector{T}, proj; σ=0.03) where {T <: Real}
    spk_filter(spk, proj; kernel=k_gaussian(σ))
end