module SpikeTriggered

import Statistics: mean

include("PSTH/PSTH.jl") # ==> binned spikes
include("Stimulus/Stimulus.jl") # ==> stimulus matrix
include("OnlineSpikes/OnlineSpikes.jl") # ==> simple event detection and removal

@doc raw"""
    STA(X, spk::Array{T}; n=10) where {T <: Real, N} -> Vector{T}

get the spike triggered average. (t0 at left)
"""
function STA(X, y::Array{T}; n=10) where {T <: Real}
    ȳ = mean(y, dims=2)
    denom = sum(ȳ)
    (N, m) = size(X)
    output = zeros(m, n)
    for tidx in 1:n
        @inbounds output[:, tidx] .= view(ȳ' * circshift(X, tidx - 1) ./ denom, :)
    end

    output[:]
end

# #TODO: STC
# function STC(X::Array{Tx, 2}, spks::Vector{Ts}) where {Tx, Ts}
#     _sta = STA(X, spks)
#     ss = X .- _sta
#     C_hat = ss * transpose(ss) ./ (size(X, 2) - 1);
#     return (eigvals(C_hat), eigvecs(C_hat))
# end

end # module
