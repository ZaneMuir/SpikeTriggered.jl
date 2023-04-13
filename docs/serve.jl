import Pkg
Pkg.activate(@__DIR__)

import LiveServer as LS
cd("build")
LS.serve()