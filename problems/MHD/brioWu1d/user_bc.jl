#---------------------------------------------------------------------------------
# Dirichlet boundary conditions of the Brio-Wu shock tube: both ends are
# held at their initial states. The left end is genuinely undisturbed at
# t = 0.2 (the left fast rarefaction's head is at x ≈ 0.13). The right end
# is not: the right-going fast wave travels at the fast speed of the right
# state, 3.75, and leaves the domain at t ≈ 0.13, after which the gas at
# x = 1 moves at u ≈ −0.24 while the boundary node is pinned at u = 0 —
# visible as a one-node glitch and a spike of the DynSGS coefficient at
# x = 1 in the output, nowhere else (the reference solution of
# reference_hll.dat uses an outflow boundary, hence its −0.24 at x = 1).
# Leaving the right end free (no value prescribed) is not an outflow
# condition for the CG discretization: the natural boundary drained the
# domain (u → −1.2 at t = 0.2, measured), so the pinned end is kept, as in
# the paper's finite-element setting.
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
