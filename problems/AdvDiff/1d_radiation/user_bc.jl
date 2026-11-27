#-------------------------------------------------------------------------------------------------
# Left: Dirichlet value, set in user_bc_dirichlet!.
# Right: radiation, -ν ∂u/∂n = εσ (u⁴ - u∞⁴). As an inward flux ν ∂u/∂n this is εσ (u∞⁴ - u⁴),
# returned by user_bc_neumann! (needs :bdy_fluxes => true). A nonlinear flux is no different
# from a linear one: it is evaluated with the current wall value at every RHS evaluation.
#-------------------------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, qe, ::TOTAL)
    if tag == "left"
        qbdy[1] = RAD_ULEFT[]
    end
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, qe, ::PERT)
    if tag == "left"
        qbdy[1] = RAD_ULEFT[] - qe[1]
    end
end

function user_bc_neumann!(F, q, q1, qe, qe1, tag, coords, τ_f, wθ, SOL_VARS_TYPE, PhysConst;
                          t = 0.0, nx = 0.0, ny = 0.0, nz = 0.0, kwargs...)
    if tag == "right"
        u = SOL_VARS_TYPE == PERT() ? q[1] + qe[1] : q[1]
        F[1] = RAD_COEFF[] * (RAD_UINF[]^4 - u^4)
    end
    return nothing
end
