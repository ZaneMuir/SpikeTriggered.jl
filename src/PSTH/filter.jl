k_gaussian(σ::Real=1) = (x) -> exp(x^2 / (- 2 * σ ^ 2)) / (σ * sqrt(2 * π))

function spk_filter(spk::Vector{T}, proj, kernel::Function) where {T <: Real}
    _psth = zeros(T, length(proj))
    for idx in 1:length(proj)
        @inbounds @fastmath _psth[idx] = sum(kernel, spk .- proj[idx])
    end
    return _psth
end

@doc raw"""
    gsl_gaussian_filter(arr::Vector{Float64}; K::Integer, α::Real=1)

The GSL gaussian kernel is defined by:
```math
G(k) = \exp\left(- \frac{k^2}{2 \sigma^2}\right)
```

for ``-(K-1)/2 \leq k \leq (K-1)/2``, and ``K`` is the size of the kernel.
The parameter ``\alpha`` specifies the number of standard deviations ``\sigma``
desired in the kernel.
```math
\sigma = \frac{K-1}{2 \alpha}
```
"""
function gsl_gaussian_filter(arr::Vector{Float64}; K::Integer, α::Real=1)
    _output = zeros(length(arr))

    _x = GSL.C.gsl_vector(length(arr), 1, pointer(arr), C_NULL, 0)
    _y = GSL.C.gsl_vector(length(_output), 1, pointer(_output), C_NULL, 0)

    _workspace = GSL.filter_gaussian_alloc(K)
    GSL.filter_gaussian(GSL.GSL_FILTER_END_PADZERO, α, 0, _x, _y, _workspace)
    GSL.filter_gaussian_free(_workspace)

    return _output
end
