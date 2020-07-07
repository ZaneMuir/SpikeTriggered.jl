module SpikeTriggered

import SDMS.get_entries
import SparseArrays: sparse, spzeros, dropzeros

include("utils.jl")
include("Stimulus/Stimulus.jl")

function STA(se::T, bspk::Vector{U}) where {T <: StimulusEnsemble, U <: Real}
    reshape(se * bspk ./ sum(bspk), se.gridSize, se.gridSize, se.temporalLength)
end

function STA(X::Array{T, 2}, bspk::Vector{U}; grid_size=16) where {T <: Real, U <: Real}
    reshape(X * bspk ./ sum(bspk), grid_size, grid_size, :)
end

end # module
