# Keep the reduced-gravity flux and wet/dry velocity treatment identical to
# ReducedGravity. DynSGS obtains g' and Hmin from this case's user_inputs.jl.
include(joinpath(@__DIR__, "..", "ReducedGravity", "user_flux.jl"))
