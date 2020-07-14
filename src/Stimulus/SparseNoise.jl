struct SparseNoise{T <: Integer} <: StimulusEnsemble
    snra::Vector{T} # s x t
    gridSize::Int64
    temporalLength::Int64
    temporalResolution::Float64
    backgroundColor::Int64 #NOTE: 0-black, 1-gray, 2-white
    multiplication::Int64 #NOTE: so far not used
    kernel::Array
    atol::Float64

    SparseNoise(snra::Vector{T}; grid_size=16, t_len=10, t_res=1/35, bgc=1, m=1, atol=1e-8) where {T <: Integer} = begin

        _k = zeros(2*(m+1-m%2), 2*(m+1-m%2))
        _k[2-m%2:end-1+m%2, 2-m%2:end-1+m%2] .= 1
        _kernel = conv(_k, ones(2, 2)/4)[2:2:end, 2:2:end]

        new{T}(snra, grid_size, t_len, t_res, bgc, m, _kernel, atol)
    end

    SparseNoise(mdb, fname::String; kwargs...) = begin

        filename = splitext(fname)[1] * ".snf"
        raw = Mongoc.open_download_stream(Mongoc.Bucket(mdb["attachments"]), filename) do mio
            read(mio)
        end
        meta = Mongoc.find_one(mdb["attachments"]["fs.files"], Mongoc.BSON("filename" => filename),
            options=Mongoc.BSON("projection"=>Dict("metadata" => 1)))["metadata"]
        dtype = Int16 #TODO: use value in meta
        SparseNoise(reinterpret(dtype, raw); kwargs...)
    end
end

function _legacy_parse_snf_array(rvec::Array{T}, grid_size::Integer) where {T}

    _aftermod = mod.(rvec .- 1, grid_size * grid_size * 2)
    _row = _aftermod .÷ (grid_size * 2) .+ 1
    _col_n = mod.(_aftermod, grid_size * 2)
    _col = _col_n .÷ 2 .+ 1
    _sign = 1 .- mod.(_col_n, 2) .* 2

    return (ids = T.(_aftermod), row = T.(_row), col = T.(_col), sign = T.(_sign))
end

Base.show(io::IO, snd::SparseNoise) = begin
    print(io, "SparseNoise{n: $(length(snd.snra)), h: $(snd.temporalLength), s: $(snd.gridSize)}")
end
# Base.eltype(snd::SparseNoise) = eltype(snd.stimulus)
Base.size(snd::SparseNoise) = (snd.gridSize * snd.gridSize * snd.temporalLength, length(snd.snra))
Base.size(snd::SparseNoise, dim::Integer) = dim > 2 ? 1 : size(snd)[dim]
Base.length(snd::SparseNoise) = length(snd.snra)

Base.getindex(snd::SparseNoise, i::Integer) = begin
    if snd.multiplication == 1
        if i < snd.temporalLength
            d = _legacy_parse_snf_array(snd.snra[1:i], snd.gridSize)

            I = d.row .+ (d.col .- 1) .* snd.gridSize .+ snd.gridSize * snd.gridSize .* ((snd.temporalLength-i):(snd.temporalLength-1))
            J = ones(eltype(snd.snra), i);
            V = d.sign
        else
            d = _legacy_parse_snf_array(snd.snra[i-snd.temporalLength+1:i], snd.gridSize)

            I = d.row .+ (d.col .- 1) .* snd.gridSize .+ snd.gridSize * snd.gridSize .* (0:(snd.temporalLength-1))
            J = ones(eltype(snd.snra), snd.temporalLength);
            V = d.sign
        end

        if snd.backgroundColor == 0
            V[V .< 0] .= 0
        elseif snd.backgroundColor == 2
            V[V .> 0] .= 0
        end

        return sparse(I, J, V, size(snd, 1), 1) |> dropzeros
    else
        dense = zeros(snd.gridSize, snd.gridSize, snd.temporalLength)

        if i < snd.temporalLength
            d = _legacy_parse_snf_array(snd.snra[1:i], snd.gridSize)
            for idx in 1:length(d.ids)
                dense[d.row[idx], d.col[idx], snd.temporalLength-i+idx] = d.sign[idx]
            end
        else
            d = _legacy_parse_snf_array(snd.snra[i-snd.temporalLength+1:i], snd.gridSize)
            for idx in 1:length(d.ids)
                dense[d.row[idx], d.col[idx], idx] = d.sign[idx]
            end
        end

        if snd.backgroundColor == 0
            dense[dense .< 0] .= 0
        elseif snd.backgroundColor == 2
            dense[dense .> 0] .= 0
        end

        _offset = div(snd.multiplication, 2)
        _tmp = conv(dense, snd.kernel)[1+_offset:end-_offset, 1+_offset:end-_offset, :]
        _tmp[abs.(_tmp) .< snd.atol] .= 0

        return sparse(_tmp[:])
    end
end

Base.getindex(snd::SparseNoise, I) = begin
    M = zeros(eltype(snd.snra), snd.temporalLength, length(I))
    for (iidx, idx) in enumerate(I)
        if idx < snd.temporalLength
            M[1:end-idx, iidx] = snd.snra[end+idx-snd.temporalLength+1:end]
            M[end-idx+1:end, iidx] = snd.snra[1:idx]
        else
            M[:, iidx] = snd.snra[idx-snd.temporalLength+1:idx]
        end
    end

    if snd.multiplication == 1
        d = _legacy_parse_snf_array(M, snd.gridSize)
        sI = d.row .+ (d.col .- 1) .* snd.gridSize .+ snd.gridSize * snd.gridSize .* (0:(snd.temporalLength-1))
        sJ = ones(size(M)) .* (1:length(I))'
        sV = d.sign

        if snd.backgroundColor == 0
            sV[sV .< 0] .= 0
        elseif snd.backgroundColor == 2
            sV[sV .> 0] .= 0
        end

        return sparse(sI[:], sJ[:], sV[:], size(snd, 1), length(I)) |> dropzeros
    else
        dense = zeros(snd.gridSize, snd.gridSize, snd.temporalLength, length(I))
        d = _legacy_parse_snf_array(M, snd.gridSize)

        if snd.backgroundColor == 0
            d.sign[d.sign .< 0] .= 0
        elseif snd.backgroundColor == 2
            d.sign[d.sign .> 0] .= 0
        end

        for idx in 1:size(d.ids, 2)
            for jdx in 1:size(d.ids, 1)
                dense[d.row[jdx, idx], d.col[jdx, idx], jdx, idx] = d.sign[jdx, idx]
            end
        end

        _offset = div(snd.multiplication, 2)
        dense = convn(dense, reshape(snd.kernel, size(snd.kernel, 1), size(snd.kernel, 2), 1, 1))[1+_offset:end-_offset, 1+_offset:end-_offset, :, :]
        dense[abs.(dense) .< snd.atol] .= 0

        return sparse(reshape(dense, :, length(I)))
    end
end

# Base.setindex!(snd::SparseNoise, v, i::Integer)
# Base.setindex!(snd::SparseNoise, v, I::Varargin{Integer, N}) where {N}

Base.firstindex(snd::SparseNoise) = 1
Base.lastindex(snd::SparseNoise) = length(snd)

Base.:*(snd::SparseNoise, bspk::Vector{T}) where {T} = snd[1:end] * bspk

# iteration
# Base.iterate(snd::SparseNoise) = size(snd, 2) > 0 ? (snd[1], 2) : nothing
# Base.iterate(snd::SparseNoise, state) = size(snd, 2) >= state ? (snd[state], state+1) : nothing
