#-------------------------------------------------------------------------------------------------
# In 1D the two boundaries are tagged "left" and "right".
#
# Dirichlet (right): set the value in user_bc_dirichlet!. The kernel imposes it on the state and
# freezes that node's tendency.
#
# Neumann (left): return the inward diffusive flux F = ν ∂u/∂n from user_bc_neumann! (needs
# :bdy_fluxes => true), and set NOTHING for that tag in user_bc_dirichlet!. Here the flux comes
# from the gradient of the exact solution projected on the outward normal nx (-1 at the left
# end), and changes in time: both t and nx are keywords of the hook.
#-------------------------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, qe, ::TOTAL)
    if tag == "right"
        qbdy[1] = n1_exact(coords[1], t)
    end
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, qe, ::PERT)
    if tag == "right"
        qbdy[1] = n1_exact(coords[1], t) - qe[1]
    end
end

function user_bc_neumann!(F, q, q1, qe, qe1, tag, coords, τ_f, wθ, SOL_VARS_TYPE, PhysConst;
                          t = 0.0, nx = 0.0, ny = 0.0, nz = 0.0, kwargs...)
    if tag == "left"
        F[1] = N1_NU[] * n1_dudx(coords[1], t) * nx      # ν ∇u · n
    end
    return nothing
end
