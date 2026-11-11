#---------------------------------------------------------------------------------
# Boundary conditions of the quarter cylinder (paper Fig. 5.1a), read by the
# quasi-static solver at every node of every boundary edge:
#
#     user_solid_bc(tag, x, y, t, nx, ny) -> NamedTuple
#
# with the edge's gmsh tag and the OUTWARD unit normal (nx, ny) of the domain.
# Fields: `ux`/`uy` prescribe that displacement component, `traction = (tx, ty)`
# prescribes σ·n; an empty NamedTuple is a traction-free surface.
#
#   "inner"   r = a    internal pressure p: σ·n = −p n. The outward normal of
#                      the solid points INTO the hole here, so −p n pushes the
#                      wall outwards, as a pressure must.
#   "outer"   r = b    traction-free (the case-bonded variant of the problem
#                      would put a displacement condition here instead).
#   "xaxis"   y = 0    symmetry: u_y = 0, u_x free.
#   "yaxis"   x = 0    symmetry: u_x = 0, u_y free.
#
# The two symmetry rollers also remove the rigid-body motions — both
# translations and the rotation — so the problem is well posed with no extra
# pinning. At the corner nodes A = (a,0) and B = (b,0) the roller and the
# pressure/free surface meet: the roller fixes u_y, the pressure acts on u_x.
#---------------------------------------------------------------------------------
function user_solid_bc(tag::String, x, y, t, nx, ny)

    if tag == "inner"
        p = cylinder_pressure(t)
        return (; traction = (-p*nx, -p*ny))
    elseif tag == "xaxis"
        return (; uy = 0.0)
    elseif tag == "yaxis"
        return (; ux = 0.0)
    elseif tag == "outer"
        return (;)
    end

    error(" # ERROR viscoelastic_cylinder2d/user_bc.jl: unexpected boundary tag \"", tag,
          "\" — the mesh should carry inner, outer, xaxis and yaxis.")
end
