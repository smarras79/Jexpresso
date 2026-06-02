"""
    Dirichlet boundary conditions for the 1D CompEuler sound-pulse case (case1).

    The left/right boundaries are held at the unperturbed Euler state
    (ρ, u, p) = (1, 0, 1).  Tags "left" and "right" are emitted by the
    1D mesh generator in src/kernel/boundaryconditions/BCs.jl.
"""
function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::TOTAL)
    PhysConst = PhysicalConst{Float64}()
    γ = PhysConst.γ

    ρ = 1.0
    u = 0.0
    p = 1.0

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

    return T(1.0), T(0.0), T(1.0)/(γ - T(1.0))
end
