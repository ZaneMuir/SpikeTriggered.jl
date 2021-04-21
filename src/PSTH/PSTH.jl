# only offers simple PSTH operations
# may consider to be an independent package if more functions are included
module PSTH

import GSL

include("histogram.jl")
include("filter.jl")

@doc """
convert PSTH to timestamps.

## Arguments
- `psth::Vector{T <: Integer}`

## Keyword Arguments
- `width`: define the bin width [default: 1]
- `sporadic`: if true, the generation would be random; otherwise, evenly divided [default: false]
- `upper`: cutoff of the PSTH bin value [default: Inf]

## Return
- `Vector{Float64}`
"""
function psth2timestamps(psth::Vector{T}; width=1, sporadic=false, upper=Inf) where {T <: Integer}
    N = sum(psth)
    output = Vector{Float64}()
    for (idx, item) in enumerate(psth)
        item == 0 && continue
        item = item < upper ? item : upper
        
        steps = if sporadic
            s = rand(item+1)
            cumsum(s ./ sum(s))[1:end-1] .* width
        else
            collect(range(width/(item+1), step=width/(item+1), length=item))
        end
        
        append!(output, steps .+ width * (idx - 1))
    end
    output
end
end
