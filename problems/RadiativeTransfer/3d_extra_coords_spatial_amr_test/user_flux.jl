function user_flux!(value::Vector, x::Vector, solL::Vector, solR::Vector, normal::Vector, iface::Int, SD::NSD_3D)
    # Upwind flux for RT problem
    # Simple test flux
    value .= 0.0
    return nothing
end
