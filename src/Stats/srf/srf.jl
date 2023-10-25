module SpatialRF

struct GaussianEllipse
    axis_major::Cdouble
    axis_minor::Cdouble
    center_x::Cdouble
    center_y::Cdouble
    rotation::Cdouble
    amplitude::Cdouble
end

area(e::GaussianEllipse) = e.axis_major * e.axis_minor * π

include("ellipse_overlap.jl")
include("schiller_overlap.jl")

function __init__()
    __init_eeover__();
end

end