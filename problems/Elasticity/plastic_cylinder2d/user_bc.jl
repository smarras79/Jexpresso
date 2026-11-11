#---------------------------------------------------------------------------------
# Boundary conditions — user_solid_bc(tag, x, y, t, nx, ny), as in
# viscoelastic_cylinder2d, except that the bore is DISPLACED, not pressed:
#
#   "inner"   r = a    u = δ(t) (x, y)/a — a radial displacement δ, no slip
#   "outer"   r = b    traction-free
#   "xaxis"   y = 0    symmetry: u_y = 0
#   "yaxis"   x = 0    symmetry: u_x = 0
#
# At the two corners of the bore the normal component of u is zero by symmetry
# AND by the bore condition. It is left to the roller alone, so that every
# constrained DOF belongs to exactly one tag: the reaction on "inner" is then
# the force the bore actually carries, and p = R_inner,x / a exactly (the
# quarter bore's resultant of a pressure p is p a in each direction).
#---------------------------------------------------------------------------------
function user_solid_bc(tag::String, x, y, t, nx, ny)

    if tag == "inner"
        P  = plastic_cylinder_properties()
        δa = bore_displacement(t)/P.a
        abs(y) < 1.0e-10 && return (; ux = δa*x)            # corner on the x axis
        abs(x) < 1.0e-10 && return (; uy = δa*y)            # corner on the y axis
        return (; ux = δa*x, uy = δa*y)
    elseif tag == "xaxis"
        return (; uy = 0.0)
    elseif tag == "yaxis"
        return (; ux = 0.0)
    elseif tag == "outer"
        return (;)
    end

    error(" # ERROR plastic_cylinder2d/user_bc.jl: unexpected boundary tag \"", tag,
          "\" — the mesh should carry inner, outer, xaxis and yaxis.")
end
