# only offers simple PSTH operations
# may consider to be an independent package if more functions are included
module PSTH

import GSL

include("histogram.jl")
include("filter.jl")
end
