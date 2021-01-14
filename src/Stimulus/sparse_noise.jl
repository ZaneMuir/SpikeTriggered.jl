
struct SparseNoise{T}
    ids::Vector{Int16}
    row::Vector{Int16}
    col::Vector{Int16}
    sign::Vector{T}
    gridsize::Integer

    SparseNoise(snra::Vector{T}, gridsize::Integer) where {T <: Integer} = begin
        _snra = parse_snra(snra, gridsize)
        new{eltype(_snra.sign)}(_snra.ids, _snra.row, _snra.col, _snra.sign, gridsize)
    end
end

function parse_snra(snra::Array{T}, grid_size::Integer) where {T <: Integer}
    _aftermod = mod.(snra .- 1, grid_size * grid_size * 2)
    _row = _aftermod .÷ (grid_size * 2) .+ 1
    _col_n = mod.(_aftermod, grid_size * 2)
    _col = _col_n .÷ 2 .+ 1
    _sign = 1 .- mod.(_col_n, 2) .* 2

    return (ids = T.(_aftermod), row = T.(_row), col = T.(_col), sign = T.(_sign))
end

Base.collect(sn::SparseNoise{T}) where {T} = begin
    N = length(sn.ids)
    output = zeros(T, N, sn.gridsize * sn.gridsize)
    for idx in 1:N
        # @inbounds output[idx, (sn.row[idx]-1)*sn.gridsize+sn.col[idx]] = sn.sign[idx]
        @inbounds output[idx, (sn.col[idx]-1)*sn.gridsize+sn.row[idx]] = sn.sign[idx]
    end
    return output
end
