#---------------------------------------------------------------------------------
# Boundary conditions of the quarter plate (Artioli et al. Fig. 5.3a):
# "displacement boundary restraints are prescribed for normal components on
# symmetry boundaries and on top and lateral boundaries".
#
#   "top"      y = H         u_y = δ(t)    the load; u_x free
#   "right"    x = L         u_x = 0       the lateral restraint; u_y free
#   "left"     x = 0, y > R  u_x = 0       symmetry
#   "bottom"   y = 0, x > R  u_y = 0       symmetry
#   "hole"     r = R         traction-free
#
# Every constrained DOF belongs to one tag only (the corners join a u_x and a
# u_y condition), so the reaction on "top" — R_top_y in solid_history.csv, or
# the nodal field ry summed over the edge — is the force the top edge carries,
# per mm of thickness (plane strain). The paper's "reaction sum" (Fig. 5.4) is
# 10 times it: their elastic slopes are exactly 10x ours while Table 1 agrees,
# i.e. a 10 mm thick plate (tools/artioli2017_figures.jl).
#---------------------------------------------------------------------------------
function user_solid_bc(tag::String, x, y, t, nx, ny)

    if tag == "top"
        return (; uy = top_displacement(t))
    elseif tag == "right" || tag == "left"
        return (; ux = 0.0)
    elseif tag == "bottom"
        return (; uy = 0.0)
    elseif tag == "hole"
        return (;)
    end

    error(" # ERROR perforated_plate2d/user_bc.jl: unexpected boundary tag \"", tag,
          "\" — the mesh should carry hole, bottom, left, right and top.")
end
