@doc raw"""
    histogram(u_arr::Vector{T}, edges::Vector{U}; kwargs ...) where {T, U}

Count histogram by providing the edges. If provides `n+1` edges, it would return `n` length histogram.
Edges would work as: left bound <= value < right bound

check out the `gsl_histogram`, which use GSL library and performed much faster.

# Arguments
- `u_arr`: a vector of values
- `edges`: a vector of edges

# Optional Keword Arguments
- `asSparse`: flag for returning result in sparse vector [default: false]
- `sorted`: flag for whether edges is sorted. If false, will use built-in `sort` first [default: false]

---
### Benchmark Test:

`histogram(rand(100_000), 0:0.01:1)`:
- julia mean time: 6.136 ms
- numpy mean time: 6.21 ms
- gsl mean time:   569.501 μs
"""
function histogram(u_arr::Vector{T}, edges::U; asSparse::Bool = false, sorted = false) where {T, U}
    if asSparse
        _h = spzeros(Int64, length(edges)-1)
    else
        _h = zeros(T, length(edges)-1)
    end

    if !sorted
        arr = sort(u_arr)
    else
        arr = u_arr
    end

    _ridx = 1
    _r_lower = edges[1]
    _r_upper = edges[2]

    for i in arr
        if _r_lower <= i < _r_upper
            @inbounds _h[_ridx] += 1
        elseif i < _r_lower
            continue
        else
            while _ridx <= length(edges)-1
                if i < edges[_ridx+1]
                    @inbounds _h[_ridx] += 1
                    break
                else
                    _ridx += 1
                end
            end
        end
    end

    return _h
end

@doc raw"""
    gsl_histogram(u_arr::Vector{T}, edges::Vector{U}; kwargs ...) where {T, U}

Count histogram by providing the edges. If provides `n+1` edges, it would return `n` length histogram.
Edges would work as: left bound <= value < right bound

# Arguments
- `u_arr`: a vector of values
- `edges`: a vector of edges

# Optional Keword Arguments
- `asSparse`: flag for returning result in sparse vector [default: false]
- `sorted`: flag for whether edges is sorted. If false, will use built-in `sort` first [default: false]

---
### Benchmark Test:

`histogram(rand(100_000), 0:0.01:1)`:
- julia mean time: 6.136 ms
- numpy mean time: 6.21 ms
- gsl mean time:   569.501 μs
"""
function gsl_histogram(u_arr::Vector{T}, edges::Vector{T}; asSparse::Bool=false, sorted::Bool=false) where {T <: Cdouble}
    !sorted && sort!(edges)

    n = length(edges)-1
    gsl_hist = GSL.histogram_alloc(n)
    GSL.histogram_set_ranges(gsl_hist, edges, n+1)

    for idx in 1:length(u_arr)
        @inbounds GSL.histogram_increment(gsl_hist, u_arr[idx])
    end

    myhist = zeros(Int, n)
    for idx in 1:n
        @inbounds myhist[idx] = GSL.histogram_get(gsl_hist, idx-1)
    end

    GSL.histogram_free(gsl_hist)

    asSparse ? sparse(myhist) : myhist
end