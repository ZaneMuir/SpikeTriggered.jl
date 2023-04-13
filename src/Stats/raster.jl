@doc """
    raster(spk, markers; head, duration, tail) -> Vector{Vector{T}}

Create rasters form spike train and onset time of stimulus.

## Arguments:
- `spk::AbstractVector{T}`: spike train
- `markers::AbstractVector`: onset time of stimulus

## Keyword Arguments:
- `head`: time included before the event markers
- `duration`: event length
- `tail`: time included after the event duration
- `offset`: flag to offset all the spike time to the onset time. [default: `true``]

## Return
- rasters as `Vector{Vector{T}}`
"""
function spike_raster(spk::AbstractVector{T}, markers::AbstractVector;
    head::Real, duration::Real, tail::Real, offset::Bool=true) where {T <: Real}

    output = Vector{T}[]
    for item in markers
        _candidates = spk[(item-head) .<= spk .< (item+duration+tail)]
        if offset
            push!(output, _candidates .- item)
        else
            push!(output, _candidates)
        end
    end
    output
end
