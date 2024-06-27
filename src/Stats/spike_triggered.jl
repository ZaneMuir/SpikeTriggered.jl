@doc raw"""
    make_strf(val; gridsize) -> Array{T, 3} SxSxT

convert strf vector into 3d spatiotemporal matrix.

When plotting with heatmap, REMEMBER to use
the transversed matrix and reversed y axis.
"""
function make_strf(val::AbstractArray; gridsize::Integer,
        _gridwidth=nothing, _gridheight=nothing)
    collect(reshape(val, gridsize, gridsize, :))
end

@doc raw"""
    hstack_strf(val; gridsize) -> Array{T, 2} Sx(SxT)

convert strf 3d matrix into horizontal matrix.

When plotting with heatmap, REMEMBER to use
the transversed matrix and reversed y axis.
"""
function hstack_strf(val::AbstractArray; kwargs...)
    _strf = make_strf(val; kwargs...)
    hstack_strf(_strf)
end

function hstack_strf(val::AbstractArray{T, 3}) where {T}
    reduce(hcat, eachslice(val; dims=3))
end

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
function spike_triggered_average(X::AbstractMatrix, y::AbstractArray; n=10, norm=true)
    T = eltype(X)
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

function _stimulus_nonlinear(video, bias)
    _video = copy(video)
    if isa(bias, Symbol)
        if bias == :on
            _video[_video .== -1] .= 0
        elseif bias == :off
            _video[_video .== 1] .= 0
        else # e.g. :diff
            nothing
        end
    else
        _video[_video .== bias] .= 0
    end
    _video
end

@doc raw"""
    spike_triggered_average_suite(stimulus, psth; kwargs...)

same inputs as for `spike_triggered_average`,
    but return a tuple of STAs from difference map, on map and off map.
"""
function spike_triggered_average_suite(stimulus, psth; scale=true, kwargs...)
    _scaler = scale ? sum(psth) : 1
    rez = (;
        diff = spike_triggered_average(stimulus, psth; norm=false, kwargs...),
        on = spike_triggered_average(
            _stimulus_nonlinear(stimulus, :on), psth; norm=false, kwargs...),
        off = spike_triggered_average(
            _stimulus_nonlinear(stimulus, :off), psth; norm=false, kwargs...),
    )
    map(x->x .* _scaler, rez)
end

function spike_triggered_average_zscore_floop(X, y; n=10, bootstrap=-1, kwargs...)
    _sta = spike_triggered_average(X, y; n, norm=false, kwargs...)
    bootstrap < 0 && (return _sta) ## simple STA

    _bootstrap_offset = if bootstrap == 0  ## full range bootstrapping
        range(n+1, step=1, length=size(X, 1)-n)
    else ## specificed bootstrap number
        _random_range = randperm(size(X, 1)-2*n)
        _random_range[1:min(length(_random_range), bootstrap)]
    end

    _bsta = zeros(eltype(X), length(_sta), length(_bootstrap_offset))
    @floop for (idx, offset) in enumerate(_bootstrap_offset)
        _bsta[:, idx] .= spike_triggered_average(circshift(X, -offset), y; n, norm=false, kwargs...)
    end
    _bsta
    _μ = mean(_bsta; dims=2)[:]
    _σ = std(_bsta; dims=2)[:]
    _zscore = (_sta .- _μ) ./ _σ
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
