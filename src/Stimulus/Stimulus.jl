abstract type StimulusEnsemble end

include("SparseNoise.jl")
include("GaussianNoise.jl")

function equal_division(ids::Vector{T}, size::Integer; total=512) where {T <: Integer}
    _train_idx = zeros(Int, length(ids) - total * size)
    _test_idx = zeros(Int, total * size)

    n_train = 0
    n_test = 0

    for idx in 1:total
        _candidates = findall(ids .== idx-1)
        rand_queue = randperm(length(_candidates))

        _train_idx[n_train+1:n_train+length(rand_queue)-size] = _candidates[rand_queue[size+1:end]]
        _test_idx[n_test+1:n_test+size] = _candidates[rand_queue[1:size]]

        n_train += length(rand_queue) - size
        n_test += size

    end

    return (_train_idx, _test_idx)