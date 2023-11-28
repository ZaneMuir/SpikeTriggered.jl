if isnothing(get(ENV, "CI", nothing))
   push!(LOAD_PATH, "../src/")

   import Pkg
   Pkg.activate(joinpath(@__DIR__, ".."))
   Pkg.resolve()
   Pkg.instantiate()
end

using Documenter, SpikeTriggered

makedocs(
         modules = [SpikeTriggered],
         sitename="SpikeTriggered.jl",
         checkdocs=:exports,
         remotes=nothing,
         format = Documenter.HTML(;
                                  repolink="/dir?ci=SpikeTriggered.jl&name=SpikeTriggered.jl",
                                  disable_git=true,
                                  edit_link=nothing,
                                 ),
         pages = [
               "Home" => "index.md",
               "Raster and PSTH" => "psth.md",
               "Spike Statistics" => "stats.md",
               "Reverse Correlation" => "reverse_correlation.md",
               "Forward Correlation" => "forward_correlation.md",
               "Index" => "func_list.md"
         ],
        )
