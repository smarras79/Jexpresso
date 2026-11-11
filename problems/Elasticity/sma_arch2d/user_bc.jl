#---------------------------------------------------------------------------------
# Boundary conditions of the arch (Artioli et al. Fig. 5.5a):
#
#   "clamp"   θ = π: y = 0, −4.5 ≤ x ≤ −3.5    u_x = u_y = 0
#   "load"    θ = 0: y = 0,  3.5 ≤ x ≤ 4.5     traction (q(t), 0): along +x,
#                                              tangential to the end face
#   "inner"   r = 3.5                          traction-free
#   "outer"   r = 4.5                          traction-free
#
# The reaction on "clamp" (R_clamp_x in solid_history.csv) balances the load:
# R_clamp_x = −q · 1 mm at every converged increment.
#
# Verification mode (JEXPRESSO_SMA_VERIFY=elastic): the clamp prescribes the
# exact curved-bar displacement and the loaded end its exact shear traction
# (user_analytic.jl), so that the elastic SEM solution must equal the exact one.
#---------------------------------------------------------------------------------
function user_solid_bc(tag::String, x, y, t, nx, ny)

    if tag == "clamp"
        if sma_verify_mode()
            ex = curved_bar_exact(x, y)
            return (; ux = ex.ux, uy = ex.uy)
        end
        return (; ux = 0.0, uy = 0.0)
    elseif tag == "load"
        return sma_verify_mode() ? (; traction = (curved_bar_end_traction(x), 0.0)) :
                                   (; traction = (sma_arch_load(t)/sma_arch_properties().thick, 0.0))
    elseif tag == "inner" || tag == "outer"
        return (;)
    end

    error(" # ERROR sma_arch2d/user_bc.jl: unexpected boundary tag \"", tag,
          "\" — the mesh should carry clamp, load, inner and outer.")
end
