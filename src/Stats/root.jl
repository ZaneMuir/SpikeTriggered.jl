# Spike Train Statistics
module Stats

import GSL
import FFTW
using Random: randperm
using Statistics: mean, var, std
using SpecialFunctions: erfc, erfcinv
using Distributions: cdf, Normal
using FLoops: @floop
using FHist: Hist1D, bincounts
using ToeplitzMatrices: Circulant

const SpikeRaster{T} = Vector{Vector{T}}
const AbstractSpikeTrain{T} = AbstractVector{T}
const AbstractPSTH{T} = AbstractVecOrMat{T}    # [nTimepoints] or [nTimepoints x nRepeats]
const AbstractStimulus{T} = AbstractMatrix{T}  # [nDimensions x nTimepoints]
const AbstractMarker{T} = AbstractVecOrMat{T}  # [nTimepoints] or [nTimepoints x nRepeats]

export SpikeRaster, AbstractSpikeTrain
export make_strf, hstack_strf
export spike_raster
export spike_histogram, spike_histogram_smoothed
export spike_triggered_average, spike_triggered_average_zscore, spike_triggered_average_suite
export get_footprint_map, get_footprint_mask
export burst_detect, burst_detect_lgn, burst_detect_lgn, burst_interpolate, split_tonic_burst, split_tonic_cardinal
export spike_train_spectrum_power

include("misc.jl")
include("raster.jl")
include("psth.jl")
include("spike_triggered_average.jl")
include("footprint.jl")
include("burst.jl")
include("spectrum.jl")

#TODO: review the functions below and make documentations
include("jpsth.jl")
include("reliability.jl")
include("correlogram.jl")
include("spike_triggered_covariance.jl")
include("srf/srf.jl")
# include("spike_triggered.jl")

end
