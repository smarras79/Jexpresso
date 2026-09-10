#---------------------------------------------------------------------------------
# Dirichlet boundary conditions of the Brio-Wu shock tube: both ends are
# held at their initial states. At the final time t = 0.1 the outermost
# waves are the left fast rarefaction's head at x ≈ 0.32 and the right fast
# rarefaction's head at x ≈ 0.87 (both move at the fast speed of their
# undisturbed state, 1.8 and 3.75), so both ends are still undisturbed.
# Running past t ≈ 0.13 lets the right wave leave the domain, after which
# the pinned right end shows a one-node glitch; a free (natural) end is not
# an outflow condition for the CG discretization and drained the domain
# when tried, so the pinned ends stay, as in the paper's finite-element
# setting.
#---------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::TOTAL)
    if tag == "left"
        ρ, u, p, By = 1.0,   0.0, 1.0,  1.0
    else  # "right"
        ρ, u, p, By = 0.125, 0.0, 0.1, -1.0
    end
    Bx = 0.75
    qbdy[1] = ρ
    qbdy[2] = ρ*u
    qbdy[3] = 0.0
    qbdy[4] = p/(γ_mhd - 1.0) + 0.5*ρ*u*u + 0.5*(Bx*Bx + By*By)
    qbdy[5] = 0.0
    qbdy[6] = Bx
    qbdy[7] = By
    qbdy[8] = 0.0
    return qbdy
end

function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::PERT)
    nothing
end

function user_bc_neumann(q::AbstractArray, gradq, coords, t, inputs)
    flux = zeros(size(q, 2), 1)
    return flux
end
