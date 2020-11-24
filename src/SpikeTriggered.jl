module SpikeTriggered

import Statistics: mean
# import SDMS.get_entries
# import Mongoc
# import SparseArrays: sparse, spzeros, dropzeros
# import LinearAlgebra: eigvals, eigvecs
# import DSP: conv
# import GSL
# import Random.randperm
# include(joinpath(@__DIR__, "../deps/FastConv/FastConv.jl"))
# convn = FastConv.convn

include("PSTH/PSTH.jl") # ==> binned spikes
include("Stimulus/Stimulus.jl") # ==> stimulus matrix
#REVIEW: include("EnsembleStimulus/Stimulus.jl") # ==> stimulus ensemble objects

@doc raw"""
    STA(X, spk::Array{T}; n=10) where {T <: Real, N} -> Vector{T}

get the spike triggered average. (t0 at left)
"""
function STA(X, y::Array{T}; n=10) where {T <: Real}
    ȳ = mean(y, dims=2)
    denom = sum(ȳ)
    (N, m) = size(X)
    output = zeros(m, n)
    for tidx in 1:n
        @inbounds output[:, tidx] .= view(ȳ' * circshift(X, tidx - 1) ./ denom, :)
    end

    output[:]
end

#REVIEW:
# function STA_ensemble(X, spk::Array{T, N}) where {T <: Real, N}
#     if N == 1
#         return X * spk ./ sum(spk)
#     elseif N == 2
#         (xsize, tlen) = size(X)
#         output = zeros(xsize)
#         @inbounds for idx in 1:size(spk, 2)
#             output += X * spk[1:tlen, idx]
#         end
#         return output ./ sum(spk)
#     end
# end

# #TODO: STC
# function STC(X::Array{Tx, 2}, spks::Vector{Ts}) where {Tx, Ts}
#     _sta = STA(X, spks)
#     ss = X .- _sta
#     C_hat = ss * transpose(ss) ./ (size(X, 2) - 1);
#     return (eigvals(C_hat), eigvecs(C_hat))
# end

end # module
