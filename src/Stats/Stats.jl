# Spike Train Statistics
module Stats

import GSL
using Statistics: mean, var, std
import ..PSTH: histogram_gsl

include("raster.jl")
include("spike_triggered.jl")
include("spectrum.jl")
include("burst.jl")
include("correlogram.jl")
include("stats.jl")

end