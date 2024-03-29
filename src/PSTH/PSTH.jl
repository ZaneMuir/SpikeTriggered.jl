# only offers simple PSTH operations
# may consider to be an independent package if more functions are included
module PSTH

import GSL
import ..Stats

include("histogram.jl")
include("filter.jl")  #TODO: review filter.jl

# @doc """
#     raster(spk::Vector{T}, markers::Vector{T}; head, duration, tail) -> Vector{Vector{T}}

# Create a raster data for given spike train and event markers (for the event start only)

# ## Arguments:
# - `spk::Vector{T}`
# - `markers::Vector{T}`

# ## Keyword Arguments:
# - head: time included before the event markers [default: 0.5]
# - duration: event length [default=1.0]
# - tail: time included after the event duration [default: 0.5]
# - norm: flag to offset all the spike time to the markers. [default: true]

# ## Return
# - `Vector{Vector{T}}`
# """
# function raster(spk::Vector{T}, markers::Vector{T}; head=0.5, duration=1, tail=0.5, norm=true) where {T <: Real} # -> Vector{Vector{T}}
#     output = Vector{Vector{T}}()
#     for item in markers
#         _roi = (item-head) .< spk .< (item+duration+tail)
#         _candidates = spk[_roi]
#         if norm
#             push!(output, _candidates .- item)
#         else
#             push!(output, _candidates)
#         end
#     end
#     output
# end

@deprecate raster(spk::Vector{T}, markers::Vector{T}; head=0.5, duration=1, tail=0.5, norm=true) where {T <: Real} Stats.spike_raster(spk, markers; head, duration, tail, offset=norm)

end