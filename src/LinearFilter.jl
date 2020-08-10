
function kernel_gaussian(H::T, sigma::U) where {T <: Integer, U <: Real}
    G = (x) -> exp(-x^2/(2 * sigma ^ 2)) / (sigma * sqrt(2 * pi))
    kernel = G.(-H:H)
end

@doc raw"""
    filter_gaussian(bpsth::Array{T, N}; H=15, sigma=5) where {T, N <= 2}

The Gaussian filter convolves the input signal with a Gaussian kernel or window.

The kernel function looks like:

```math
G(k) = \frac{\exp(-\frac{1}{2} \cdot (\frac{k}{\sigma})^2)}{\sigma \sqrt{2 \pi}}
```
where ``-H \le k \le H``. It would make the window size as ``2 H + 1``

"""
function filter_gaussian(bpsth::Array{T, N}; H=15, sigma=5) where {T, N}
    kernel = kernel_gaussian(H, sigma)
    output = zeros(T, size(bpsth))

    if sigma == 0
        return bpsth
    end

    if N == 1
        output[:] = conv(bpsth, kernel)[H+1:end-H][:] #TODO: use fastconv.convn(::Array{T, N}, ::Array{T, N})
    elseif N == 2
        output[:] = conv(bpsth, kernel)[H+1:end-H, :][:]
    else
        throw(ArgumentError("bpsth must be 1d or 2d array"))
    end
    return output
end