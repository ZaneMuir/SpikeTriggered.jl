
@doc raw"""
    get_footprint_map(zscore::AbstractArray; alpha=0.01, collapse=false, gridsize=nothing)

create a footprint map based on the q-values
using Benjamini–Hochberg procedure.

color codes as:
- `1`: ON response
- `0`: NO response
- `-1`: OFF response
""" 
function get_footprint_map(zscore::AbstractArray; alpha=0.01,
        collapse=false, gridsize=nothing)
    @assert !(collapse && (gridsize isa Nothing)) "please specify gridsize to get collapsed footprint map"
    
    _mask = get_footprint_mask(zscore; alpha)
    _map = map(x->Int(sign(x)), _mask .* zscore)

    if collapse
        map(x->begin
                _idx = findfirst(abs.(x) .> 0)
                _idx isa Nothing ? 0 : Int(x[_idx])
                end, eachslice(make_strf(_map; gridsize); dims=(1, 2)))
    else
        _map
    end
end

@doc raw"""
    get_footprint_mask(zscore::AbstractArray; alpha=0.01)

create a binary footprint mask based on the q-values
using Benjamini–Hochberg procedure.
"""
function get_footprint_mask(zscore::AbstractArray; alpha=0.01)
    get_footprint_mask_from_qvalue(zscore_fdr_qvalue(zscore); alpha)
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
function zscore_fdr_qvalue(zscore::AbstractArray;)

    Pz = zscore_pvalue(zscore) .* 2  # 2 tails
    M = length(Pz)
    k = sortperm(Pz[:])
    qval = Pz .* M ./ reshape(k, size(Pz)...)

    (qval, k)
end

@doc raw"""
    get_footprint_mask_from_qvalue(rez; alpha=0.01)

create a binary footprint mask based on the q-values.
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

## # https://en.wikipedia.org/wiki/False_discovery_rate
## function _zscore_FDR_threshold(zscore::AbstractArray{T}, α::AbstractVector) where {T}
##     # threshold from false discovery rate test
##     z_cdf = normcdf.(- abs.(zscore)) .* 2
##     volumn = length(z_cdf)
##     z_cdf_norm = sort(z_cdf) ./ (1:volumn) .* volumn
##     # sort(z_cdf)
## 
##     _rez = fill(T(Inf), length(α))
##     for (idx, val) in enumerate(α)
##         _thresh_idx = findlast(z_cdf_norm .< val)
##         isnothing(_thresh_idx) && continue
## 
##         _thresh_cdf = z_cdf_norm[_thresh_idx] / volumn * _thresh_idx
##         _thresh_z = abs(norminv(_thresh_cdf/2))
##         _rez[idx] = _thresh_z
##     end
##     _rez
## end
## _zscore_FDR_threshold(zscore::AbstractArray, α::T) where {T <: Real} = zscore_FDR_threshold(zscore, T[α])[1]
## 
## function _calc_footprint(zval, thresh::T; val=0x01, flag=0x80) where {T <: Real}
##     rez = zeros(UInt8, length(zval))
##     _val = val .| flag
##     rez[abs.(zval) .> thresh] .= _val
##     rez
## end
## 
## @doc raw"""
##     spike_triggered_footprint(stimulus, psth; sta_n=10, bootstrap_n=100, α=())
## 
## Footprint as a general approach to quantifying RF size, based on statistical methods used for fMRI.
## With bootstrapping dataset, use false discovery rate test to identify significant spatiotemporal pixels.
## Check out Soto-Sanchez et al 2017 for more details.
## 
## """
## function spike_triggered_footprint(stimulus, psth; sta_n=10, bootstrap_n=100, α::NTuple{aN, T}=()) where {aN, T}
##     @assert 0 <= aN <= 3 "α should have length between 0 to 3."
## 
##     ## 1 get STA
##     _sta = spike_triggered_average_suite(stimulus, psth; n=sta_n)
## 
##     ## 2 bootstrapping
##     N = min(length(psth), bootstrap_n)
##     _sta_perm = map(
##         x->spike_triggered_average_suite(stimulus, circshift(psth, x); n=sta_n),
##         randperm(length(psth))[1:N]
##     )
## 
##     ## 3 sta zscore
##     sta_z_diff = let
##         _sta_shuffle = stack(map(x->x.diff, _sta_perm))
##         _mu = mean(_sta_shuffle; dims=2)[:]
##         _sig = std(_sta_shuffle; dims=2)[:]
##         (_sta.diff .- _mu) ./ _sig
##     end
##     sta_z_on = let
##         _sta_shuffle = stack(map(x->x.on, _sta_perm))
##         _mu = mean(_sta_shuffle; dims=2)[:]
##         _sig = std(_sta_shuffle; dims=2)[:]
##         (_sta.on .- _mu) ./ _sig
##     end
##     sta_z_off = let
##         _sta_shuffle = stack(map(x->x.off, _sta_perm))
##         _mu = mean(_sta_shuffle; dims=2)[:]
##         _sig = std(_sta_shuffle; dims=2)[:]
##         (_sta.off .- _mu) ./ _sig
##     end
## 
##     ## 4 significant stixel from false discovery rate test
##     A = Float64[0.05, 0.01, 0.001, α...]
##     get_full_footprint = (sta_zscore, val)->begin
##         _thresh = _zscore_FDR_threshold(sta_zscore, A)
##         _footprint = reduce((a, b) ->a .| b, map((a)->_calc_footprint(sta_zscore, a[2]; val, flag=0x01 << (8 - a[1])), enumerate(_thresh)); init=zeros(UInt8, length(sta_zscore)))
##         (; footprint=_footprint, threshold=_thresh, alpha=A)
##     end
## 
##     (;
##         meta = (; bootstrapN = N),
##         sta = _sta,
##         diff = get_full_footprint(sta_z_diff, 0x03),
##         on = get_full_footprint(sta_z_on, 0x01),
##         off = get_full_footprint(sta_z_off, 0x02),
##     )
## end
## 
## @doc raw"""
##     get_footprint(raw, threshIdx=1; s=16, collapse=false)
## 
## get the spatial footprint from the footprint calculation.
## """
## function get_footprint(raw, threshIdx=1; s=16, collapse=false)
##     _fp = zeros(Bool, length(raw.footprint))
##     _flag = 0x01 << (8 - threshIdx)
##     _roi = findall(x->x .& _flag == _flag, raw.footprint)
##     _fp[_roi] .= true
##     _fp = reshape(_fp, s, s, :)
## 
##     if collapse
##         sum(_fp, dims=3)[:, :] .> 0
##     else
##         _fp
##     end
## end
