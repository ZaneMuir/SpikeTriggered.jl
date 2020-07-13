char2luminance(x; dtype::Type=Float64) = dtype(x / 255 * 2 - 1)

struct GaussianNoise{Tv <: Real} <: StimulusEnsemble
    video::Array{Tv, 3}
    temporalLength::Int
    temporalResolution::Float64
    converter::Function

    GaussianNoise(video::Array{T, 3}; t_len::Int, t_res=1/35, converter::Function=char2luminance) where {T} = begin
        new{T}(video, t_len, t_res, converter)
    end
end

Base.show(io::IO, gnd::GaussianNoise) = print(io, "gaussian noise design - size: $(size(gnd.video)), t_len: $(gnd.temporalLength)")

Base.size(gnd::GaussianNoise) = (size(gnd.video, 1) * size(gnd.video, 2) * gnd.temporalLength, size(gnd.video, 3))
Base.size(gnd::GaussianNoise, dim::Integer) = dim > 2 ? 1 : size(gnd)[dim]
Base.length(gnd::GaussianNoise) = size(gnd.video, 3)

Base.getindex(gnd::GaussianNoise, i::Integer) = begin
    #     frame = zeros(size(gnd, 1))
        if i < gnd.temporalLength
            V = zeros(size(gnd, 1))
            _step = size(gnd.video, 1) * size(gnd.video, 2)
            V[1:_step*(gnd.temporalLength - i)] = gnd.video[:,:,end-gnd.temporalLength+i+1:end][:]
            V[_step*(gnd.temporalLength - i)+1:end] = gnd.video[:,:,1:i][:]
        else
            T = (i-gnd.temporalLength+1):i
            V = gnd.video[:,:,T][:]
        end
        gnd.converter.(V)
    end

Base.getindex(gnd::GaussianNoise, I) = begin
    Vs = zeros(size(gnd, 1), length(I))
    for (iidx, idx) in enumerate(I)
        Vs[:, iidx] = gnd[idx]
    end
    Vs
end

Base.firstindex(gnd::GaussianNoise) = 1
Base.lastindex(gnd::GaussianNoise) = length(gnd)

Base.:*(gnd::GaussianNoise, bspk::Vector{T}) where {T} = gnd[1:end] * bspk