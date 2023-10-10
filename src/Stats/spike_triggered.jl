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

#TODO: CR
@doc raw"""
    spike_triggered_covariance(X::AbstractMatrix{T}, y::AbstractArray{T}; n::Integer=10, verbose=true) where {T, U} -> Matrix{T}
"""
function spike_triggered_covariance(X::AbstractMatrix{T}, y::AbstractArray{T}; n::Integer=10, verbose=true) where {T}
    _Nx = size(X, 1)
    
    _sta = spike_triggered_average(X, y; n, norm=false)
    _sta_c = _sta ./ sum(abs2, _sta)
    
    _stc = zeros(T, length(_sta), length(_sta))
    
    _y = mean(y, dims=2)[:]
    
    pbar = PM.Progress(length(_y))
    
    for idx in 1:length(_y)
        _ind = (((idx):(-1):(idx-n+1)) .- 1) .% _Nx .+ 1
        _ind[_ind .< 1] .= _ind[_ind .< 1] .+ _Nx
        
        _s = X[_ind, :][:]
        _s = _s .- _s' * _sta .* _sta_c

        _stc .+= (_s * _s') .* y[idx]
        
        verbose && PM.next!(pbar)
    end
    _stc ./= sum(_y) - 1
    _stc
end