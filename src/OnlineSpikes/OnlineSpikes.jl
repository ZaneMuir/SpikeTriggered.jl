module OnlineSpikes

"""
    detect_cross_pnt(arr; thresh, gap=1, dir=:up, ref=nothing)
detect data rise or data down points.
# Arguments
- `arr`: a vector of traces.
- `thresh`: threshold
"""
function detect_cross_pnt(arr::Vector{T}; thresh::T, gap::Int=1, dir::Symbol=:up, ref=nothing, kwargs...) where {T}
    above_thresh_idx = findall(arr .> thresh)
    above_gap_idx_idx = findall(diff(above_thresh_idx) .> gap)
    
    if dir == :up
        candidate_idx = above_thresh_idx[above_gap_idx_idx .+ 1]
        if above_thresh_idx[1] > 1
            candidate_idx = [above_thresh_idx[1]; candidate_idx]
        end
        
    elseif dir == :down
        candidate_idx = above_thresh_idx[above_gap_idx_idx]
        if above_thresh_idx[end] < length(arr)
            candidate_idx = [candidate_idx; above_thresh_idx[end]]
        end
    else
        throw(ArgumentError("unknown direction: $(dir)"))
    end
    
    if isnothing(ref)
        return candidate_idx
    else
        return ref[candidate_idx]
    end
end

function detect_events(trace; kwargs...)
    cross_pnt = detect_cross_pnt(trace; kwargs...)
    detect_events(trace, cross_pnt; kwargs...)
end

function detect_events(trace, cross_pnt::Vector{T}; method::String, win_size, kwargs...) where {T}
    criteria = if method == "maxima"
        findmax
    elseif method == "minima"
        findmin
    end
    
    output = zeros(Int64, length(cross_pnt))
    for (idx, item) in enumerate(cross_pnt)
        _, offset = criteria(trace[item-win_size:item+win_size])
        output[idx] = item - win_size + offset
    end
    
    output
end

include("event_removal.jl")

end