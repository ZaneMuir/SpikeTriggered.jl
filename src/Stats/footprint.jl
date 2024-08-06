
@doc raw"""
    get_footprint_map([func_pval,] raw::AbstractArray; alpha=0.01, fdr_c=1, collapse=false, gridsize=nothing)

Create a footprint map based on the q-values from `func_pval(raw)`
using Benjamini–Hochberg procedure; `alpha` as the threshold, fdr_c as the constant.

```math
P_k \leq \frac{k}{N \cdot C} \alpha
```

For arbiturary join distribution of p-values,
one could set `fdr_c` as `ln(N)+1/N+γ (Euler constant)`.
(using `benjamini_hochberg_constant(N)`)

If `collapse` is `true` and `gridsize` is specified, a 2d time-collapsed
footprint map will be returned.

color codes as:
- `1`: ON response
- `0`: NO response
- `-1`: OFF response
"""
function get_footprint_map(func_pval::Function, raw::AbstractArray; collapse=false, gridsize=nothing, kwargs...)
    @assert !(collapse && (gridsize isa Nothing)) "please specify gridsize to get collapsed footprint map"

    _mask = get_footprint_mask(func_pval, raw; kwargs...)
    _map = sign.(_mask .* raw)

    if collapse
        #NOTE: polarity is set as the first sign of footprint
        map(x->begin
                _idx = findfirst(abs.(x) .> 0)
                _idx isa Nothing ? 0 : x[_idx]
                end, eachslice(make_strf(_map; gridsize); dims=(1, 2)))
    else
        _map
    end
end
get_footprint_map(raw::AbstractArray; kwargs...) = get_footprint_map(zscore_pvalue, raw; kwargs...)

@doc raw"""
    get_footprint_mask([func_pval::Function,] raw::AbstractArray; alpha=0.01, fdr_c=1)

Create a binary footprint mask based with Benjamini–Hochberg procedure using arbiturary q-values
converted from `raw` statistic values by `func_pval`.

In most cases, one should use `get_footprint_map` which can also generate 3d mask
with addtional polarity information.

If `func_pval` is not set, `zscore_pvalue` will be used, which assumes
inputs are zscore from standard normal distribution.

`alpha` if the threshold for false discovery rate, default as `0.01`.

`fdr_c` is the constant value for correcting Benjamini–Hochberg procedure with joint distributions among pvalues.
"""
function get_footprint_mask(func_pval::Function, raw::AbstractArray; alpha=0.01, fdr_c=1)
    pval = func_pval(raw)
    qval, qk = benjamini_hochberg_qvalue(pval; C=fdr_c)
    get_footprint_mask_from_qvalue(qval, invperm(qk[:]); alpha) #NOTE: invperm only takes a vector.
end
get_footprint_mask(raw::AbstractArray; kwargs...) = get_footprint_mask(zscore_pvalue, raw; kwargs...)

#TODO: docs
#NOTE: when N > 50, the difference is less than 0.01
#NOTE: when N > 2000, the `sum` starts to be significantly slower than the approximation
benjamini_hochberg_constant(N::Integer) = N>50 ? log(N) + Base.MathConstants.eulergamma + 1/N : sum(x->1/x, 1:N)

@doc raw"""
    benjamini_hochberg_qvalue(pvalue; C=1) -> (q::Array, k::Array)

False discovery rate test for any list of pvalues.
This function will return the `q` value and the corresponding `k` index.

```math
P_k \leq \frac{k}{N} \alpha
```

```math
q_k = \frac{N}{k} P_k \leq \alpha
```

`k` is generated from `invperm(sortperm(...))`.
"""
function benjamini_hochberg_qvalue(pvalue::AbstractArray; C=1)
    N = length(pvalue)
    k = reshape(invperm(sortperm(pvalue[:])), size(pvalue)...)
    qval = pvalue .* (N * C) ./ k

    (qval, k)
end

@doc raw"""
    zscore_pvalue(zscore::AbstractArray; two_tailed=true)

get pvalue from standard normal distribution
"""
function zscore_pvalue(zscore::AbstractArray; two_tailed=true)
    p = cdf(Normal(), -1 .* abs.(zscore))
    two_tailed ? 2 .* p : p
end

@doc raw"""
    get_footprint_mask_from_qvalue(qval, k; alpha=0.01)

create a binary footprint mask based on the q-values and the corresponding permutation vector.
Using Benjamini–Hochberg procedure.
"""
function get_footprint_mask_from_qvalue(qval, k; alpha=0.01)
    _cutoff = findlast(qval[k] .<= alpha)

    _mask = zeros(Bool, size(qval))
    if !isnothing(_cutoff)
        _mask[k[1:_cutoff]] .= true
    end

    _mask
end
