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

function interp_burst(burst_iti::Vector{T}; n=13, iterp_t=nothing) where {T <: Real}
    N = length(burst_iti)
    iti_list = Cdouble.(burst_iti)
    output = zeros(Cdouble, n)

    xa = collect(range(0, stop=1, length=N))
    iterp_type = if iterp_t == :linear
        GSL.gsl_interp_linear
    elseif iterp_t == :polynomial
        GSL.gsl_interp_polynomial
    elseif iterp_t == :cspline
        GSL.gsl_interp_cspline
    else
        GSL.gsl_interp_steffen
    end

    #TODO: maybe use GSL.gsl_spline instead.
    interp_obj = GSL.interp_alloc(iterp_type, N)
    accel_obj = GSL.interp_accel_alloc()

    GSL.interp_init(interp_obj, xa, iti_list, N)
    
    for (idx, i) in enumerate(range(0, stop=1, length=n))
        output[idx] = GSL.interp_eval(interp_obj, xa, iti_list, Cdouble(i), accel_obj)
    end

    GSL.interp_accel_free(accel_obj)
    GSL.interp_free(interp_obj)
    
    output
end