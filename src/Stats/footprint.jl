
@doc raw"""
    get_footprint_map(zscore::AbstractArray; alpha=0.01, fdr_c=1, collapse=false, gridsize=nothing)

Create a footprint map based on the q-values from z-scores
using Benjamini–Hochberg procedure and `alpha` as the threshold.

```math
P_k \leq \frac{k}{N} \alpha
```

For arbiturary join distribution of p-values,
one could set `fdr_c` as `ln(N)+1/N+γ (Euler constant)`.

If `collapse` is `true` and `gridsize` is specified, a 2d time-collapsed
footprint map will be returned.

color codes as:
- `1`: ON response
- `0`: NO response
- `-1`: OFF response
"""
function get_footprint_map(zscore::AbstractArray; alpha=0.01,
        collapse=false, gridsize=nothing, fdr_c=1)
    @assert !(collapse && (gridsize isa Nothing)) "please specify gridsize to get collapsed footprint map"

    _mask = get_footprint_mask(zscore; alpha, fdr_c)
    _map = map(x->Int(sign(x)), _mask .* zscore)

    if collapse
        #NOTE: polarity is set as the first sign of footprint
        map(x->begin
                _idx = findfirst(abs.(x) .> 0)
                _idx isa Nothing ? 0 : Int(x[_idx])
                end, eachslice(make_strf(_map; gridsize); dims=(1, 2)))
    else
        _map
    end
end

@doc raw"""
    get_footprint_mask(zscore::AbstractArray; alpha=0.01, fdr_c=1)

Create a binary footprint mask based on the q-values
using Benjamini–Hochberg procedure.

In most cases, user should use `get_footprint_map` which can also generate 3d mask
with addtional polarity information.
"""
function get_footprint_mask(zscore::AbstractArray; alpha=0.01, fdr_c=1)
    get_footprint_mask_from_qvalue(zscore_fdr_qvalue(zscore; C=fdr_c); alpha)
end

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
    zscore_fdr_qvalue(zscore::AbstractArray) -> (qval::AbstractVector, k::Vector{Int})

Computed the corresponding q-values of Z-scores for each spatiotemporal voxel,
using Benjamini–Hochberg multiple-hypothesis testing.

Assuming all voxels are independent.

Use `get_footprint_from_qvalue` to generate footprints.

Returns:
- `qval`: the q-values
- `k`: the permutation vector of the zscore, based on the 2-tailed p-values of normal distribution.
"""
function zscore_fdr_qvalue(zscore::AbstractArray; C=1)

    Pz = zscore_pvalue(zscore) .* 2  # 2 tails
    (qval, k) = benjamini_hochberg_qvalue(Pz; C)
    return (qval, invperm(k[:]))
end

@doc raw"""
    get_footprint_mask_from_qvalue(qval, k; alpha=0.01)

create a binary footprint mask based on the q-values and the corresponding permutation vector.
Using Benjamini–Hochberg procedure.
"""
function get_footprint_mask_from_qvalue(rez::Tuple; kwargs...)
    get_footprint_mask_from_qvalue(rez...; kwargs...)
end

function get_footprint_mask_from_qvalue(qval, k; alpha=0.01)
    _cutoff = findlast(qval[k] .<= alpha)

    _mask = zeros(Bool, size(qval))
    if !isnothing(_cutoff)
        _mask[k[1:_cutoff]] .= true
    end

    _mask
end
