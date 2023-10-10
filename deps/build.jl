
println("building EEOver")
run(Cmd(`make julia-lib`; dir=joinpath(@__DIR__, "eeover-c")))