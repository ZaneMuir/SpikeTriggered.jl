@doc raw"""
    spike_triggered_average(X::AbstractMatrix{T}, y::AbstractArray{T}; n=10, norm=true) where {T} -> Vector{T}

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
function spike_triggered_average(X::AbstractMatrix{T}, y::AbstractArray{T}; n=10, norm=true) where {T}
    ȳ = mean(y, dims=2)
    denom = sum(abs, ȳ)
    (N, m) = size(X)
    output = zeros(T, m, n)
    for tidx in 1:n
        @inbounds output[:, tidx] .= view(ȳ' * circshift(X, tidx - 1) ./ denom, :)
    end

    if norm
        output[:] ./ maximum(abs, output)
    else
        output[:]
    end
end
