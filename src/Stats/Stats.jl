# Spike Train Statistics
module Stats

import GSL
using Statistics: mean, var, std

include("raster.jl")
include("psth.jl")
include("jpsth.jl")
include("spike_triggered.jl")
include("spectrum.jl")
include("burst.jl")
include("correlogram.jl")
include("stats.jl")

end