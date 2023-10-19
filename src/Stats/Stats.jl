# Spike Train Statistics
module Stats

import GSL
import FFTW
using Random: randperm
using Statistics: mean, var, std
using SpecialFunctions: erfc, erfcinv

include("raster.jl")
include("psth.jl")
include("jpsth.jl")
include("spike_triggered.jl")
include("spectrum.jl")
include("burst.jl")
include("correlogram.jl")
include("stats.jl")
include("ellipse_overlap.jl")
include("footprint.jl")

end