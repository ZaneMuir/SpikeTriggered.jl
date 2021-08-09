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

@doc """
    raster(spk::Vector{T}, markers::Vector{T}; head, duration, tail) -> Vector{Vector{T}}

Create a raster data for given spike train and event markers (for the event start only)

## Arguments:
- `spk::Vector{T}`
- `markers::Vector{T}`

## Keyword Arguments:
- head: time included before the event markers [default=0.5]
- duration: event length [default=1.0]
- tail: time included after the event duration [default=0.5]

## Return
- `Vector{Vector{T}}`
"""
function raster(spk::Vector{T}, markers::Vector{T}; head=0.5, duration=1, tail=0.5) where {T <: Real} # -> Vector{Vector{T}}
    output = Vector{Vector{T}}()
    for item in markers
        _roi = (item-head) .< spk .< (item+duration+tail)
        _candidates = spk[_roi]
        push!(output, _candidates .- item)
    end
    output
end
