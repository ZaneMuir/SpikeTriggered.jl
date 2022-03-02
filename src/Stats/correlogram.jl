function correlogram(target::Vector{T}, reference::Vector{T}; bin_size=1e-2, window_size=300, is_auto=false) where {T}
    _K_half = div(window_size, 2)
    _bin_boundary = (-_K_half:_K_half) .* bin_size
    output = zeros(Int64, window_size)
    for pt in reference
        _offset = target .- pt
        _offset = is_auto ? _offset[_offset .!= 0] : _offset
        output .+= SpikeTriggered.PSTH.histogram_gsl(_offset, _bin_boundary)
    end
    _bin_boundary[2:end] .- bin_size/2, output
end

#TODO:
# function correlogram()

autocorrelogram(ts_vec; kwargs...) = correlogram(ts_vec, ts_vec; is_auto=true, kwargs...)
crosscorrelogram(args...; kwargs...) = correlogram(args...; is_auto=false, kwargs...)