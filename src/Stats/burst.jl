@doc raw"""
    detect_burst(spk::Vector{T}; t_silence=0.07, t_isi=0.03, nofs=nothing) where {T <: Real} -> Vector{Vector{T}}

Burst firing identification based on the rules:

- a burst would followed/preceded by a silence lasting `t_silence` time
- each spike spaced at most `t_isi` time, otherwise terminated
- (optional) at least `nofs[1]` spikes within the first `nofs[2]` time

References:
- [Vaingankar et al 2012](https://doi.org/10.3389/fnint.2012.00118)

## Arguments
- `spk::Vector{T}`: spike train vector. one should make sure it is orderd.

## Keyword Arguments:
- `t_silence`
- `t_isi`
- `nofs`

## Returns
- `burst_list::Vector{Vector{T}}`
"""
function detect_burst(spk::Vector{T}; t_silence=0.07, t_isi=0.03, nofs::Union{Nothing, Tuple{Integer, Real}}=nothing) where {T <: Real}
    
    # initialization
    burst_list = Vector{Vector{T}}()
    current_burst = Vector{T}()
    flag_in_burst = false

    # given the spike numbers are usually large, it is okay to skip the first and the last bursts, if any.
    for idx in 2:(length(spk)-1)
        if flag_in_burst
            
            # continue
            if spk[idx+1] - spk[idx] <= t_isi
                push!(current_burst, spk[idx])
            
            # the last one
            else
                push!(current_burst, spk[idx])
                flag_in_burst = false
                push!(burst_list, current_burst)
            end
            
        # seek for potential burst start
        else
            if spk[idx] - spk[idx-1] >= t_silence && spk[idx+1] - spk[idx] <= t_isi
                flag_in_burst = true
                current_burst = Vector{T}([spk[idx]])
            end
        end
    end
    
    if !isnothing(nofs)
        output = Vector{Vector{T}}()
        while !isempty(burst_list)
            candidate = popfirst!(burst_list)
            if length(candidate) >= nofs[1]
                if sum(diff(candidate[1:nofs[1]])) <= nofs[2]
                    push!(output, candidate)
                end
            end
        end
        output
    else
        burst_list
    end
end