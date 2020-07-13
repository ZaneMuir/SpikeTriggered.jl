module SpikeTriggered

import SDMS.get_entries
import SparseArrays: sparse, spzeros, dropzeros
import DSP: conv
include(joinpath(@__DIR__, "../deps/FastConv/FastConv.jl"))
convn = FastConv.convn

include("utils.jl")
include("Stimulus/Stimulus.jl")

function STA(se::T, bspk::Vector{U}; flatten=false) where {T <: StimulusEnsemble, U <: Real}
    if flatten
        reshape(se * bspk ./ sum(bspk), se.gridSize, se.gridSize, se.temporalLength)
    else
        se * bspk ./ sum(bspk)
    end
end

function STA(X::Array{T, 2}, bspk::Vector{U}; grid_size=16, flatten=false) where {T <: Real, U <: Real}
    if flatten
        reshape(X * bspk ./ sum(bspk), grid_size, grid_size, :)
    else
        X * bspk ./ sum(bspk)
    end
end

function STA(X, bspks::Array{T, 2}; kwargs...) where {T <: Real}
    n = size(bspks, 2)
    result = zeros(size(X, 1), n)

    for idx in 1:n
        result[:, idx] = STA(X, bspks[:, idx]; flatten=true, kwargs...)[:] .* sum(bspks[:, idx])
    end

    sum(result, dims=2)[:] ./ sum(bspks)
end

end # module
