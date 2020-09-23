char2luminance(x; dtype::Type=Float64) = dtype(x / 255 * 2 - 1)

struct GaussianNoise{Tv <: Real} <: StimulusEnsemble
    video::Array{Tv, 3}
    temporalLength::Int
    temporalResolution::Float64
    converter::Function

    GaussianNoise(video::Array{T, 3}; t_len::Int=10, t_res=1/35, converter::Function=char2luminance) where {T} = begin
        new{T}(video, t_len, t_res, converter)
    end

    GaussianNoise(mdb, fname; kwargs...) = begin
        raw = Mongoc.open_download_stream(Mongoc.Bucket(mdb["attachments"]), fname) do mio
            read(mio)
        end
        meta = Mongoc.find_one(mdb["attachments"]["fs.files"], Mongoc.BSON("filename" => fname))["metadata"]

        GaussianNoise(reshape(raw, Int64(meta["width"]), Int64(meta["height"]), Int64(meta["length"])); kwargs...)
    end
end

# rendering the GN video as dense matrix
Base.collect(gnd::GaussianNoise; tlen::Union{Nothing, Integer}=nothing, repeats::Integer=1) = begin
    if isnothing(tlen)
        video = zeros(size(gnd) .* (1, repeats))
        tlen = size(gnd, 2)
    else
        video = zeros(size(gnd, 1), tlen * repeats)
    end

    _tmp = gnd[1:tlen]
    for idx in 1:repeats
        video[:, (idx-1)*tlen+1:tlen*idx] .= _tmp
    end

    video
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

# deprecated
Base.:*(gnd::GaussianNoise, bspk::Vector{T}) where {T} = gnd[1:end] * bspk