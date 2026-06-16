"""
    Dirichlet boundary conditions for the 1D Sod shock-tube case (sod1d).

    Until the rarefaction reaches the left wall and the shock reaches the
    right wall (well past t = 0.2 for the canonical case), the boundary
    states are simply the constant left/right initial conditions.
"""
function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::TOTAL)
    PhysConst = PhysicalConst{Float64}()
    γ = PhysConst.γ

    if tag == "left"
        ρ, u, p = 1.000, 0.0, 1.0
    else  # "right"
        ρ, u, p = 0.125, 0.0, 0.1
    end

    qbdy[1] = ρ
    qbdy[2] = ρ*u
    qbdy[3] = p/(γ - 1.0) + 0.5*ρ*u*u
    return qbdy
end

function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::PERT)
    nothing
end

function user_bc_neumann(q::AbstractArray, gradq, coords, t, inputs)
    flux = zeros(size(q, 2), 1)
    return flux
end

function user_bc_dirichlet_gpu(q, qe, coords, t, lpert)
    T = eltype(q)
    PhysConst = PhysicalConst{Float64}()
    γ = T(PhysConst.γ)
    # Left state by default; the CPU dispatch above handles tag-based selection.
    return T(1.0), T(0.0), T(1.0)/(γ - T(1.0))
end
