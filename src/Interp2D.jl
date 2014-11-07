module Interp2D

const lib624 = joinpath(dirname(@__FILE__), "../deps/src/lib624.so")
if (dlopen_e(lib624) == C_NULL)
    error("Interp2D not properly installed. Run Pkg.build(\"Interp2D\")")
end

include("backend.jl")
include("frontend.jl")

export Linear2DInterpolator, Cubic2DInterpolator, evaluate, evaluate!

end