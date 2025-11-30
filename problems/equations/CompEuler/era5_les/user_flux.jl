"""
    ERA5 LES - Flux Functions

    Standard compressible Euler fluxes for LES.

    Author: Jexpresso Development Team
    Date: 2025-11-30
"""

function user_flux!(F::SubArray, G::SubArray, H::SubArray,
                    SD::NSD_3D,
                    q::AbstractArray,
                    qe::AbstractArray,
                    mesh::St_mesh,
                    ::CL, ::TOTAL, ::CompEuler;
                    neqs=7)
    PhysConst = PhysicalConst{Float64}()
    γ_gas = PhysConst.γ_dry

    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]
    ρθ = q[5]
    ρqt = q[6]
    ρql = q[7]
    P = q[end]

    u = ρu/ρ
    v = ρv/ρ
    w = ρw/ρ

    # Fluxes for compressible Euler equations
    F[1] = ρu
    F[2] = ρu*u + P
    F[3] = ρu*v
    F[4] = ρu*w
    F[5] = ρu*ρθ/ρ
    F[6] = ρu*ρqt/ρ
    F[7] = ρu*ρql/ρ

    G[1] = ρv
    G[2] = ρv*u
    G[3] = ρv*v + P
    G[4] = ρv*w
    G[5] = ρv*ρθ/ρ
    G[6] = ρv*ρqt/ρ
    G[7] = ρv*ρql/ρ

    H[1] = ρw
    H[2] = ρw*u
    H[3] = ρw*v
    H[4] = ρw*w + P
    H[5] = ρw*ρθ/ρ
    H[6] = ρw*ρqt/ρ
    H[7] = ρw*ρql/ρ
end

function user_flux!(F::SubArray, G::SubArray, H::SubArray,
                    SD::NSD_3D,
                    q::AbstractArray,
                    qe::AbstractArray,
                    mesh::St_mesh,
                    ::CL, ::PERT, ::CompEuler;
                    neqs=7)
    PhysConst = PhysicalConst{Float64}()

    # Perturbation formulation
    ρ  = q[1] + qe[1]
    ρu = q[2] + qe[2]
    ρv = q[3] + qe[3]
    ρw = q[4] + qe[4]
    ρθ = q[5] + qe[5]
    ρqt = q[6] + qe[6]
    ρql = q[7] + qe[7]
    P = q[end]

    u = ρu/ρ
    v = ρv/ρ
    w = ρw/ρ

    F[1] = ρu
    F[2] = ρu*u + P
    F[3] = ρu*v
    F[4] = ρu*w
    F[5] = ρu*ρθ/ρ
    F[6] = ρu*ρqt/ρ
    F[7] = ρu*ρql/ρ

    G[1] = ρv
    G[2] = ρv*u
    G[3] = ρv*v + P
    G[4] = ρv*w
    G[5] = ρv*ρθ/ρ
    G[6] = ρv*ρqt/ρ
    G[7] = ρv*ρql/ρ

    H[1] = ρw
    H[2] = ρw*u
    H[3] = ρw*v
    H[4] = ρw*w + P
    H[5] = ρw*ρθ/ρ
    H[6] = ρw*ρqt/ρ
    H[7] = ρw*ρql/ρ
end
