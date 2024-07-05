# Statistics

@doc raw"""
    dispersion_index(psth; corrected=true, kwargs...) -> Real

Also known as the `Fano factor`.

```math
D = \frac{\sigma^2}{\mu}
```

## optional arguments:
- `dims`: calculate the index over dimension `dims`.

Check out more on [wikipedia](https://en.wikipedia.org/wiki/Index_of_dispersion).
"""
function dispersion_index(psth; corrected=true, kwargs...)
    var(psth; corrected, kwargs...) ./ mean(psth; kwargs...)
end

# gaussian distribution cumulative density, single sided, use `cdf(Normal(), x)` instead
normcdf(x::T) where {T <: Real} = erfc(- x / sqrt(T(2))) / 2

# inverse of normcdf: use `quantile(Normal(), x)` instead
norminv(x::T) where {T <: Real} = erfcinv(2 * x) * (- sqrt(T(2)))

@doc raw"""
    make_strf(val; gridsize) -> Array{T, 3} SxSxT

convert strf vector into 3d spatiotemporal matrix.

When plotting with heatmap, REMEMBER to use
the transversed matrix and reversed y axis.
"""
function make_strf(val::AbstractArray; gridsize::Integer,
        _gridwidth=nothing, _gridheight=nothing)
    #TODO: check which axis is width and which is height
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