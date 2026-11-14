"""
    conserve_shallow_water_volume!(RHS, params)

Opt-in (`:swe_conserve_volume => true`): remove the basin mean of the assembled
continuity tendency, so that the total volume ∫H dA = Σᵢ Mᵢ Hᵢ is constant to
round-off whatever the continuity equation's source does.

Called in `_build_rhs!` after the global DSS and before the division by the
lumped mass, where `RHS[i, 1]` holds ∫ψᵢ ∂H/∂t dA. Its sum over the nodes is the
total volume tendency: the flux divergence contributes only its boundary flux,
which the walls make zero, so what is left is ∫ S₁ dA — e.g. the deep-water
entrainment of `ReducedGravity_entrain`. Subtracting w̄·Mᵢ with
w̄ = Σᵢ RHS[i,1] / Σᵢ Mᵢ removes it as a uniform detrainment w̄ (m/s) over the
basin, without changing anything else in the tendency.

Shared interface nodes hold the full assembled value on every rank, so the two
sums run over owned nodes only (`gip2owner`, 0-based MPI ranks) and are reduced
in one in-place `Allreduce!` of two numbers. Nothing is allocated (checked in
`test/test_shallow_water_volume.jl`).

Assumes the lumped (diagonal) mass matrix and no Dirichlet condition on H, CPU.
Mass added outside the RHS — the `:swe_limiter` clip — is not compensated.
"""
const _SWE_VOLUME_BUF = zeros(Float64, 2)

function conserve_shallow_water_volume!(RHS::AbstractMatrix, params)
    M = params.M
    M isa AbstractVector ||
        error(":swe_conserve_volume needs the lumped (diagonal) mass matrix of the inexact LGL quadrature.")
    return _conserve_volume_kernel!(RHS, M, params.mesh.gip2owner,
                                    Int(params.mesh.npoin), get_mpi_comm())
end

# Function barrier: every argument concretely typed.
function _conserve_volume_kernel!(RHS::AbstractMatrix, M::AbstractVector,
                                  owner::AbstractVector, npoin::Int, comm)
    rank  = MPI.Comm_rank(comm)
    s_rhs = 0.0
    s_m   = 0.0
    @inbounds for ip in 1:npoin
        if owner[ip] == rank
            s_rhs += RHS[ip, 1]
            s_m   += M[ip]
        end
    end
    buf = _SWE_VOLUME_BUF
    buf[1] = s_rhs
    buf[2] = s_m
    MPI.Allreduce!(buf, MPI.SUM, comm)
    w = buf[1] / buf[2]
    @inbounds for ip in 1:npoin
        RHS[ip, 1] -= w * M[ip]
    end
    return w
end

function validate_shallow_water_volume_option(inputs, params)
    get(inputs, :swe_conserve_volume, false) == true || return nothing
    params.backend isa CPU || error(":swe_conserve_volume is implemented for the CPU backend only.")
    params.AD isa ContGal  || error(":swe_conserve_volume is implemented for continuous Galerkin only.")
    get(inputs, :ladapt, false) == true && error(":swe_conserve_volume does not support AMR (:ladapt).")
    params.M isa AbstractVector ||
        error(":swe_conserve_volume needs the lumped (diagonal) mass matrix of the inexact LGL quadrature.")
    return nothing
end
