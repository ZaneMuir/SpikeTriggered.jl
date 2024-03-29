# @doc raw"""
#     histogram(u_arr::Vector{T}, edges::Vector{U}; kwargs ...) where {T, U}

# Count histogram by providing the edges. If provides `n+1` edges, it would return `n` length histogram.
# Edges would work as: left bound <= value < right bound

# check out the `gsl_histogram`, which use GSL library and performed much faster.

# # Arguments
# - `u_arr`: a vector of values
# - `edges`: a vector of edges

# # Optional Keword Arguments
# - `sorted`: flag for whether edges is sorted. If false, will use built-in `sort` first [default: false]

# ---
# ### Benchmark Test:

# `histogram(rand(100_000), 0:0.01:1)`:
# - julia mean time: 5.957 ms
# - numpy mean time: 5.29 ms
# - gsl mean time:   467.721 μs
# """
# function histogram(u_arr::Vector{T}, edges; sorted = false) where {T}
#     if !sorted
#         arr = sort(u_arr)
#     else
#         arr = u_arr
#     end

#     _h = zeros(Int64, length(edges)-1)

#     _ridx = 1
#     _r_lower = edges[1]
#     _r_upper = edges[2]

#     for i in arr
#         if _r_lower <= i < _r_upper
#             @inbounds _h[_ridx] += 1
#         elseif i < _r_lower
#             continue
#         else
#             while _ridx <= length(edges)-1
#                 if i < edges[_ridx+1]
#                     @inbounds _h[_ridx] += 1
#                     break
#                 else
#                     _ridx += 1
#                 end
#             end
#         end
#     end

#     return _h
# end

# @doc raw"""
#     histogram_gsl(u_arr::Vector{T}, edges) where {T <: Real}

# Count histogram by providing the edges. If provides `n+1` edges, it would return `n` length histogram.
# Edges would work as: left bound <= value < right bound

# # Arguments
# - `u_arr`: a vector of values
# - `edges`: a vector of edges

# ---
# ### Benchmark Test:

# `histogram(rand(100_000), 0:0.01:1)`:
# - julia mean time: 5.957 ms
# - numpy mean time: 5.29 ms
# - gsl mean time:   467.721 μs
# """
# function histogram_gsl(u_arr::Vector{T}, edges) where {T <: Real}
#     u_arr = Cdouble.(u_arr)
#     edges = Cdouble.(edges)

#     edges = sort(edges)

#     n = length(edges)-1
#     gsl_hist = GSL.histogram_alloc(n)
#     GSL.histogram_set_ranges(gsl_hist, edges, n+1)

#     for idx in 1:length(u_arr)
#         @inbounds GSL.histogram_increment(gsl_hist, u_arr[idx])
#     end

#     myhist = zeros(Int, n)
#     for idx in 1:n
#         @inbounds myhist[idx] = GSL.histogram_get(gsl_hist, idx-1)
#     end

#     GSL.histogram_free(gsl_hist)

#     return myhist
# end

@deprecate histogram_gsl(u_arr::Vector{T}, edges) where {T <: Real} Stats.spike_histogram(u_arr, edges)

@doc """
    psth2timestamps(psth::Vector{T}; width=1, sporadic=false, upper=Inf) where {T <: Integer} -> Vector{Float64}

convert PSTH to timestamps.

## Arguments
- `psth::Vector{T <: Integer}`

## Keyword Arguments
- `width`: define the bin width [default: 1]
- `sporadic`: if true, the generation would be random (using `rand`); otherwise, evenly divided [default: false]
- `upper`: cutoff of the PSTH bin value [default: Inf]

## Return
- `Vector{Float64}`
"""
function psth2timestamps(psth::Vector{T}; width=1, sporadic=false, upper=Inf) where {T <: Integer}
    N = sum(psth)
    output = Vector{Float64}()
    for (idx, item) in enumerate(psth)
        item == 0 && continue
        item = item < upper ? item : upper

        steps = if sporadic
            s = rand(item+1)
            cumsum(s ./ sum(s))[1:end-1] .* width
        else
            collect(range(width/(item+1), step=width/(item+1), length=item))
        end

        append!(output, steps .+ width * (idx - 1))
    end
    output
end

#XXX: optimization required!
#NOTE: use spike filter instead.
@doc raw"""
    psth2trace(hist::Vector{T}, roi; σ=0.03, freq=35) where {T} -> Vector{T}

convert histgram into trace-like data, with a Gaussian filter.

## Arguments
- `hist::Vector{T}`: histogram vector
- `roi`: range or a vector of timepoints.

## Keyword Arguments:
- `σ`: sigma of the gaussian kernel (in seconds) [default: 0.03]
- `freq`: equivalent sampling frequency (Hz) of the histogram. [default: 35]

## Returns:
- `Vector{T}`: with the same length as the `roi`.
"""
function psth2trace(hist::Vector{T}, roi; σ=0.03, freq=35) where {T}
    ts = range(1/freq/2, step=1/freq, length=length(hist))
    output = zeros(length(roi))
    # k_func = SpikeTriggered.PSTH.k_gaussian(σ)
    k_func = (x) -> exp(x^2 / (-2 * σ ^ 2)) # kernel function

    Threads.@threads for idx in 1:length(output)
        @inbounds output[idx] = sum(hist .* k_func.(ts .- roi[idx] .+ 1/freq/2))
    end

    output
end