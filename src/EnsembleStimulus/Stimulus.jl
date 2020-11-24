module Ensembles

abstract type StimulusEnsemble end

include("SparseNoise.jl")
include("GaussianNoise.jl")
# include("ScrambleSN.jl")
end
