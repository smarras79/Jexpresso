function user_bc!(value::AbstractArray, x::Vector, y::Vector, z::Vector, normal::Vector, iface::Int, SD::NSD_3D)
    # Boundary condition for RT problem
    # Zero inflow/outflow for simplicity in this test case
    value .= 0.0
    return nothing
end
