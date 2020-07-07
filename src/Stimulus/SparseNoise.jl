struct SparseNoise{T <: Integer} <: StimulusEnsemble
    snra::Vector{T} # s x t
    gridSize::Int64
    temporalLength::Int64
    temporalResolution::Float64
    backgroundColor::Int64 #NOTE: 0-black, 1-gray, 2-white
    multiplication::Int64 #NOTE: so far not used

    SparseNoise(snra::Vector{T}; grid_size=16, t_len=10, t_res=1/35, bgc=1, m=1) where {T <: Integer} = begin
        new{T}(snra, grid_size, t_len, t_res, bgc, m)
    end

    SparseNoise(mdb, fname::String; kwargs...) = begin
        d = get_entries(mdb["inventory"]["attachments"], """{"filename": "$(fname)"}""", limit=1)[1]
        SparseNoise(Int64.(d["data"]["snra"]); kwargs...)
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
    if i < snd.temporalLength
        d = _legacy_parse_snf_array(snd.snra[1:i], snd.gridSize)

        I = d.row .+ (d.col .- 1) .* snd.gridSize .+ snd.gridSize * snd.gridSize .* ((10-i):9)
        J = ones(eltype(snd.snra), i);
        V = d.sign
    else
        d = _legacy_parse_snf_array(snd.snra[i-9:i], snd.gridSize)

        I = d.row .+ (d.col .- 1) .* snd.gridSize .+ snd.gridSize * snd.gridSize .* (0:9)
        J = ones(eltype(snd.snra), 10);
        V = d.sign
    end

    if snd.backgroundColor == 0
        V[V .< 0] .= 0
    elseif snd.backgroundColor == 2
        V[V .> 0] .= 0
    end

    sparse(I, J, V, size(snd, 1), 1) |> dropzeros
end

Base.getindex(snd::SparseNoise, I) = begin
    M = zeros(eltype(snd.snra), snd.temporalLength, length(I))
    for (iidx, idx) in enumerate(I)
        if idx < snd.temporalLength
            M[1:end-idx, iidx] = snd.snra[end+idx-snd.temporalLength+1:end]
            M[end-idx+1:end, iidx] = snd.snra[1:idx]
        else
            M[:, iidx] = snd.snra[idx-9:idx]
        end
    end

    d = _legacy_parse_snf_array(M, snd.gridSize)
    sI = d.row .+ (d.col .- 1) .* snd.gridSize .+ snd.gridSize * snd.gridSize .* (0:9)
    sJ = ones(size(M)) .* (1:length(I))'
    sV = d.sign

    if snd.backgroundColor == 0
        sV[sV .< 0] .= 0
    elseif snd.backgroundColor == 2
        sV[sV .> 0] .= 0
    end

    sparse(sI[:], sJ[:], sV[:], size(snd, 1), length(I)) |> dropzeros
end

# Base.setindex!(snd::SparseNoise, v, i::Integer)
# Base.setindex!(snd::SparseNoise, v, I::Varargin{Integer, N}) where {N}

Base.firstindex(snd::SparseNoise) = 1
Base.lastindex(snd::SparseNoise) = length(snd)

Base.:*(snd::SparseNoise, bspk::Vector{T}) where {T} = snd[1:end] * bspk

# iteration
# Base.iterate(snd::SparseNoise) = size(snd, 2) > 0 ? (snd[1], 2) : nothing
# Base.iterate(snd::SparseNoise, state) = size(snd, 2) >= state ? (snd[state], state+1) : nothing
