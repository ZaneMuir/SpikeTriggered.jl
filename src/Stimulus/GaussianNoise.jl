char2luminance(x; dtype::Type=Float64) = dtype(x / 255 * 2 - 1)

struct GaussianNoise{Tv <: Real} <: StimulusEnsemble
    Video::Array{Tv, 3}
    TemporalLength::Int
    TemporalResolution::Float64
    Converter::Function

    GaussianNoise(video::Array{T, 3}, t_len::Int; t_res=1/35, converter::Function=char2luminance) where {T} = begin
        new{T}(video, t_len, t_res, converter)
    end
end

Base.show(io::IO, gnd::GaussianNoise) = print(io, "gaussian noise design - size: $(size(gnd.Video)), t_len: $(gnd.TemporalLength)")

Base.size(gnd::GaussianNoise) = (size(gnd.Video, 1) * size(gnd.Video, 2) * gnd.TemporalLength, size(gnd.Video, 3))
Base.size(gnd::GaussianNoise, dim::Integer) = dim > 2 ? 1 : size(gnd)[dim]
Base.length(gnd::GaussianNoise) = size(gnd.Video, 3)

Base.getindex(gnd::GaussianNoise, i::Integer) = begin
    #     frame = zeros(size(gnd, 1))
        if i < gnd.TemporalLength
            V = zeros(size(gnd, 1))
            _step = size(gnd.Video, 1) * size(gnd.Video, 2)
            V[1:_step*(gnd.TemporalLength - i)] = gnd.Video[:,:,end-gnd.TemporalLength+i+1:end][:]
            V[_step*(gnd.TemporalLength - i)+1:end] = gnd.Video[:,:,1:i][:]
        else
            T = (i-gnd.TemporalLength+1):i
            V = gnd.Video[:,:,T][:]
        end
        gnd.Converter.(V)
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