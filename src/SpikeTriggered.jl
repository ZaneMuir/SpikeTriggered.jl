module SpikeTriggered

import SDMS.get_entries
import Mongoc
import SparseArrays: sparse, spzeros, dropzeros
import LinearAlgebra: eigvals, eigvecs
import DSP: conv
import GSL
import Random.randperm
include(joinpath(@__DIR__, "../deps/FastConv/FastConv.jl"))
convn = FastConv.convn

include("utils.jl")
include("LinearFilter.jl")
include("Stimulus/Stimulus.jl")

@doc raw"""
    STA(X, spk::Array{T, N}) where {T <: Real, N} -> Vector{T}

get the spike triggered average.
"""
function STA(X, spk::Array{T, N}) where {T <: Real, N}
    if N == 1
        return X * spk ./ sum(spk)
    elseif N == 2
        (xsize, tlen) = size(X)
        output = zeros(xsize)
        @inbounds for idx in 1:size(spk, 2)
            output += X * spk[1:tlen, idx]
        end
        return output ./ sum(spk)
    end
end

#NOTE: this function could change very dramatically
#TODO: optimize and make the results more useful
function STC(X::Array{Tx, 2}, spks::Vector{Ts}) where {Tx, Ts}
    _sta = STA(X, spks)
    ss = X .- _sta
    C_hat = ss * transpose(ss) ./ (size(X, 2) - 1);
    return (eigvals(C_hat), eigvecs(C_hat))
end

end # module
