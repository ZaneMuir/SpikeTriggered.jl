module SpikeTriggered

include("PSTH/PSTH.jl") # ==> binned spikes
include("Stimulus/Stimulus.jl") # ==> stimulus matrix
include("Waveforms/Waveforms.jl") # ==> simple event detection and removal
include("Stats/Stats.jl") # ==> spike related statistics: sta, stc, burst, power spectrum, etc.

end # module
