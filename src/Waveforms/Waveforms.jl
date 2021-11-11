module Waveforms
import GSL
using Statistics: median

include("events.jl")
include("filter.jl")

#TODO: ERP function
function trace_event_related(trace::Vector{T}, markers::Vector{T}; head=0.5, duration=1, tail=0.5, freq=10_000) where {T} # -> Vector{Vector{T}}
    output = Vector{Vector{T}}()
    _tidx = range(1/freq/2, step=1/freq, length=length(trace))
    for item in markers
        _roi = (item - head) .< _tidx .<= (item + duration + tail)
        push!(output, trace[_roi])
    end
    output
end

end