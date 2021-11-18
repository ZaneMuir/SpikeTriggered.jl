# Spike Train Statistics
module Stats

import GSL
using Statistics: mean

include("spectrum.jl")
include("burst.jl")

@doc raw"""
    spike_triggered_average(X, y::Array{T}; n=10) where {T<:Real} -> Vector

get the spike triggered average from the stimulus matrix and the spike trains.

## Arguments:
- `X`: stimulus matrix that supports `:*` and `circshift` functions, with general shape of [nTimePoints x nDimensions].
In most cases, try to use the built-in Array type.
- `y::Array{T}`: spike train as column vectors. STA will be calculated with the trial average.

## Keyword Arguments:
- `n`: how many timepoints to look back [default: 10].

## Returns:
- `Vector`: Vector of the type of `X`; flattened version of the STA matrix [nDimensions x n]. (NOTE: t0 at index `1`.)
"""
function spike_triggered_average(X, y::Array{T}; n=10, norm=true) where {T <: Real}
    ȳ = mean(y, dims=2)
    denom = sum(abs, ȳ)
    (N, m) = size(X)
    output = zeros(m, n)
    for tidx in 1:n
        @inbounds output[:, tidx] .= view(ȳ' * circshift(X, tidx - 1) ./ denom, :)
    end

    if norm
        output[:] ./ maximum(abs, output)
    else
        output[:]
    end
end

# #TODO: STC
# function STC(X::Array{Tx, 2}, spks::Vector{Ts}) where {Tx, Ts}
#     _sta = STA(X, spks)
#     ss = X .- _sta
#     C_hat = ss * transpose(ss) ./ (size(X, 2) - 1);
#     return (eigvals(C_hat), eigvecs(C_hat))
# end

end