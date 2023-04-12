push!(LOAD_PATH, "../src/")

import Pkg
Pkg.activate(@__DIR__)
Pkg.resolve()
Pkg.instantiate()
using Documenter, SpikeTriggered

makedocs(
    sitename = "SpikeTriggered.jl",
    pages = [
        "Home" => "index.md",
        "Raster and PSTH" => "psth.md",
        "Reverse Correlation" => "reverse_correlation.md",
        "Forward Correlation" => "forward_correlation.md",
        "Index" => "func_list.md"
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
)