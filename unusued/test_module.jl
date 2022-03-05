module Coordinates

export Coord, distance_from_origin
#=
The parameter `T` will be replaced by the type
used to hold the values of `x` and `y`. 
Presumably, you want the `x` and `y` of one coordinate
to share a single type (maybe Float64 or Int32).

When you construct an actual Coord with Float64 values,
this happens:

julia> point = Coord(1.0, 3.0)
Coord(1.0, 3.0)
julia> point.y
3.0
=#
struct Coord{T}
    x::T 
    y::T 
end 

function distance_from_origin(point::Coord{T}) where {T}
    return sqrt(point.x * point.x + point.y * point.y)
end

end # module Coordinates
