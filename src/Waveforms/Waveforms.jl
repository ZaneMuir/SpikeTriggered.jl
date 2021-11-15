module Waveforms
import GSL
using Statistics: median

include("events.jl")
include("filter.jl")

#TODO: ERP function
function trace_event_related(trace::Vector{T}, markers::Vector{T}; head=0.5, duration=1, tail=0.5, freq=10_000) where {T} # -> Vector{Vector{T}}
    output = Vector{Vector{T}}(undef, length(markers))
    
    _N = length(trace)
    _marker_idx = round.(Int64, markers .* freq) .+ 1
    _start_offset = round(Int64, head * freq) - 1
    _end_offset = round(Int64, (duration+tail) * freq) - 1
    
    for (idx, item) in enumerate(_marker_idx)
        _roi = max(1, item-_start_offset):min(_N, item+_end_offset)
        @inbounds output[idx] = trace[_roi]
    end

    output
end

end