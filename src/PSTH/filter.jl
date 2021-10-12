k_gaussian(σ::Real=1) = (x) -> exp(x^2 / (- 2 * σ ^ 2)) / (σ * sqrt(2 * π))

function spk_filter(spk::Vector{T}, proj, kernel::Function) where {T <: Real}
    _psth = zeros(T, length(proj))
    for idx in 1:length(proj)
        @inbounds @fastmath _psth[idx] = sum(kernel, spk .- proj[idx])
    end
    return _psth
end
