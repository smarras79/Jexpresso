#Rain
include("../friction.jl")
function user_source_friction(SD, T, q::Array, npoin::Int64)

    return friction(SD, q, 0.033, npoin; law=1)   
 
end
