module Stimulus

import Base

include("sparse_noise.jl")

function process_first_of_paired_markers(m::Vector{T}) where {T}
    @assert length(m) % 2 == 0 "markers have to be paired!"
    m[1:2:end-1]
end

#TODO: doc
function apply_video_nonlinearity(video, bias=0)
    _video = copy(video)
    _video[_video .== bias] .= 0
    _video
end

end
