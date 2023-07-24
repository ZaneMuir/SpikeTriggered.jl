#TODO: doc
function spike_jpsth_raw(x_raster::Vector{Vector{T}}, y_raster::Vector{Vector{T}}, roi::AbstractVector) where {T <: Real}
    @assert(length(x_raster) == length(y_raster), "JPSTH inputs must have same length.")
    jpsth = zeros(Int, length(roi)-1, length(roi)-1)
    for idx in eachindex(x_raster)
        x_psth = spike_histogram(x_raster[idx], roi)
        y_psth = spike_histogram(y_raster[idx], roi)
        jpsth .+= x_psth * y_psth'
    end
    jpsth
end

function spike_jpsth_shifted(x_raster::Vector{Vector{T}}, y_raster::Vector{Vector{T}}, roi::AbstractVector) where {T <: Real}
    _raw = spike_jpsth_raw(x_raster, y_raster, roi)
    shift_predictor = zeros(Int, size(_raw)..., length(x_raster))
    for idx in eachindex(x_raster)
        shift_predictor[:, :, idx] .= spike_jpsth_raw(x_raster, circshift(y_raster, idx-1), roi)
    end
    _raw .- mean(shift_predictor, dims=3)[:, :, 1], shift_predictor
end