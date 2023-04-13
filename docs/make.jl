if isnothing(get(ENV, "CI", nothing))
    push!(LOAD_PATH, "../src/")

    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.resolve()
    Pkg.instantiate()
end

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

if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "github.com/ZaneMuir/SpikeTriggered.jl.git",
    )
end