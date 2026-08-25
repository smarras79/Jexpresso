#=============================================================================
 schur_kernel.jl -- the bespoke scalar element sweeps for the Schur operator H,
                    replacing the two full five-field operator applications the
                    reference form uses.

 WHY, WITH THE MEASUREMENT THAT FORCED IT
 -----------------------------------------
 schur.jl builds H out of two calls to `hevi_apply_A!` on purpose, so that it
 cannot disagree with the operator it preconditions, and says an optimised
 kernel "can be added later and checked against THIS". This is that kernel.

 The reference form is not merely un-optimised, it is SLOWER PER ITERATION THAN
 THE FIVE-FIELD SOLVE IT REPLACES. Measured two ways:

     schur_H! / hevi_apply_A!      2.038x   (mock, p=4, 18785 points)
     matvec / f_imp                2.079x   (production, 25 ranks, 2.1M points)

 So the "scalar" matvec costs two five-field applications. On the 20x20x80
 bubble that made matvec 66% of the step and turned the whole reduction into a
 4.7% wash, even though it halves the iteration count and cuts the
 preconditioner 4.9x, the orthogonalisation 12.3x and the MPI reduce 8.9x.
 Every part of the projection held except this one, and this one is the
 dominant term.

 WHAT THE REFERENCE FORM SPENDS IT ON
 ------------------------------------
 `_hevi_A_elements_full!` computes, for every node, nine derivative chains per
 implicit field -- 45 chains -- because it does not know which of its five rows
 the caller wants. schur_grad! wants three momentum rows that are just the
 gradient of ONE scalar; schur_divW! wants one continuity row. Two applications
 is 90 chains to obtain what these sweeps get in 12:

   gradient    dP/dxi, dP/deta, dP/dzeta computed ONCE and contracted with the
               metric three ways                                    ->  3 chains
   divergence  three components, three chains each                  ->  9 chains

 plus none of the flux staging, none of the four rows that are discarded, and a
 fifth of the memory traffic through F/G/H.

 EXACTLY WHAT THESE MUST REPRODUCE, read off the reference path rather than
 re-derived (the reason the reference stays in the tree):

   schur_grad!   sets V[:,stheta] = P./beta, so beta*dRHOtheta == P, making the
                 momentum fluxes F[.,su] = G[.,sv] = H[.,sw] = P and every other
                 flux zero. Row q = su then sees only dF/dx, so
                 rhs_el[su] = -wJac * dP/dx, and gx = -W[:,su] after DSS and the
                 mass division. Net:      gx = M^-1 DSS( wJac * dP/dx ).
                 NO WALL MASKING: free slip kills the normal MASS flux, not the
                 pressure term in the normal MOMENTUM flux.

   schur_divW!   sets the momentum slots, so F[.,srho] = mx etc. with mx, my, mz
                 wall-masked, and row q = srho sees the three mass-flux
                 derivatives. Net:        out = M^-1 DSS( wJac * div(W m) ).
                 MASKING IS APPLIED, and it is the mask that makes this the same
                 divergence the continuity row takes.

 The advective Theta row does not enter either: it replaces the Theta row only,
 and neither sweep reads it.

 THE FUNCTION BARRIER IS NOT OPTIONAL. St_metrics and St_mesh carry untyped
 fields, so `metrics.dzetadz` infers as `Any` and every index inside an element
 loop becomes a boxed dispatch -- measured at 60x on the operator this file is
 speeding up. The arrays are passed positionally for that reason, exactly as
 `hevi_apply_A!` passes them to its own kernels.
=============================================================================#

"""
    SchurKernel

Scratch for the fast sweeps. Allocated once; nothing here is per-call.

`rel3`/`R3` carry the gradient's three components and `rel1`/`R1` the
divergence's one, in the same element-then-assemble layout `hevi_apply_A!` uses,
so the DSS, the halo exchange and the mass division are the SAME functions
applied to a narrower array rather than a second implementation of them.

Its own storage rather than a view into `op.rhs_el`: the operator's scratch is
in fact free while these run, but borrowing it would couple this file to the
promise that no `hevi_apply_A!` is ever in flight around it, and the whole
workspace is ~5 MB per rank against a Krylov basis of 85 MB.
"""
struct SchurKernel
    rel3::Array{Float64, 5}      # (nelem, ngl, ngl, ngl, 3)  gradient
    R3::Matrix{Float64}          # (npoin, 3)
    rel1::Array{Float64, 5}      # (nelem, ngl, ngl, ngl, 1)  divergence
    R1::Matrix{Float64}          # (npoin, 1)
    Pel::Array{Float64, 3}       # element-local scalar
    Mx::Array{Float64, 3}        # element-local masked momentum
    My::Array{Float64, 3}
    Mz::Array{Float64, 3}
    vaux::Vector{Float64}
end

function SchurKernel(nelem::Int, ngl::Int, npoin::Int)
    SchurKernel(zeros(Float64, nelem, ngl, ngl, ngl, 3),
                zeros(Float64, npoin, 3),
                zeros(Float64, nelem, ngl, ngl, ngl, 1),
                zeros(Float64, npoin, 1),
                zeros(Float64, ngl, ngl, ngl),
                zeros(Float64, ngl, ngl, ngl),
                zeros(Float64, ngl, ngl, ngl),
                zeros(Float64, ngl, ngl, ngl),
                zeros(Float64, npoin))
end

SchurKernel(params) = SchurKernel(Int(params.mesh.nelem), Int(params.mesh.ngl),
                                  Int(params.mesh.npoin))

#-----------------------------------------------------------------------------
# Element sweeps. Positional arrays: see the barrier note in the header.
#-----------------------------------------------------------------------------
"""
    _schur_grad_elements!(rel, P, Pel, conn, dpsi, om, Je, metrics..., nelem, ngl)

`rel[iel,i,j,k,1:3] = wJac * grad(P)`, the weak gradient before assembly.

THREE CHAINS, NOT NINE. dP/dxi, dP/deta and dP/dzeta do not depend on which
Cartesian component is being formed, so they are computed once per node and
contracted with the three metric rows. The reference path recomputes an
equivalent set once per implicit field.
"""
function _schur_grad_elements!(rel, P, Pel, conn, dψ, ω, Je,
                               dξdx, dξdy, dξdz,
                               dηdx, dηdy, dηdz,
                               dζdx, dζdy, dζdz,
                               nelem::Int, ngl::Int)
    @inbounds for iel = 1:nelem
        # Gather once into a contiguous element block. The chains below read it
        # three ways (m,j,k / i,m,k / i,j,m) and the strided reads are what this
        # loop is bound by -- the same finding that made branching out the zero
        # fluxes worthless in the five-field kernel.
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            Pel[i, j, k] = P[conn[iel, i, j, k]]
        end
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ωJac = ω[i] * ω[j] * ω[k] * Je[iel, i, j, k]
            dPξ = 0.0; dPη = 0.0; dPζ = 0.0
            for m = 1:ngl
                dPξ += dψ[m, i] * Pel[m, j, k]
                dPη += dψ[m, j] * Pel[i, m, k]
                dPζ += dψ[m, k] * Pel[i, j, m]
            end
            rel[iel, i, j, k, 1] = ωJac * (dPξ * dξdx[iel,i,j,k] +
                                           dPη * dηdx[iel,i,j,k] +
                                           dPζ * dζdx[iel,i,j,k])
            rel[iel, i, j, k, 2] = ωJac * (dPξ * dξdy[iel,i,j,k] +
                                           dPη * dηdy[iel,i,j,k] +
                                           dPζ * dζdy[iel,i,j,k])
            rel[iel, i, j, k, 3] = ωJac * (dPξ * dξdz[iel,i,j,k] +
                                           dPη * dηdz[iel,i,j,k] +
                                           dPζ * dζdz[iel,i,j,k])
        end
    end
    return nothing
end

"""
    _schur_div_elements!(rel, vx, vy, vz, wall, wallx, wally, Mx, My, Mz, ...)

`rel[iel,i,j,k,1] = wJac * div(W v)`, the weak divergence of the WALL-MASKED
field, before assembly.

The mask is applied at the flux, node by node, exactly where
`_hevi_A_elements_full!` applies it -- free slip means the normal MASS flux
vanishes, and it is that mask which makes this the same divergence the
continuity row takes, which is in turn what lets the Schur elimination close.
"""
function _schur_div_elements!(rel, vx, vy, vz, wall, wallx, wally,
                              Mx, My, Mz, conn, dψ, ω, Je,
                              dξdx, dξdy, dξdz,
                              dηdx, dηdy, dηdz,
                              dζdx, dζdy, dζdz,
                              nelem::Int, ngl::Int)
    @inbounds for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = conn[iel, i, j, k]
            Mx[i, j, k] = wallx[ip] ? 0.0 : vx[ip]
            My[i, j, k] = wally[ip] ? 0.0 : vy[ip]
            Mz[i, j, k] = wall[ip]  ? 0.0 : vz[ip]
        end
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ωJac = ω[i] * ω[j] * ω[k] * Je[iel, i, j, k]
            axξ = 0.0; axη = 0.0; axζ = 0.0
            ayξ = 0.0; ayη = 0.0; ayζ = 0.0
            azξ = 0.0; azη = 0.0; azζ = 0.0
            for m = 1:ngl
                dξm = dψ[m, i]; dηm = dψ[m, j]; dζm = dψ[m, k]
                axξ += dξm * Mx[m, j, k]; axη += dηm * Mx[i, m, k]; axζ += dζm * Mx[i, j, m]
                ayξ += dξm * My[m, j, k]; ayη += dηm * My[i, m, k]; ayζ += dζm * My[i, j, m]
                azξ += dξm * Mz[m, j, k]; azη += dηm * Mz[i, m, k]; azζ += dζm * Mz[i, j, m]
            end
            # Only the DIAGONAL contractions: d(mx)/dx + d(my)/dy + d(mz)/dz.
            rel[iel, i, j, k, 1] =
                ωJac * ((axξ * dξdx[iel,i,j,k] + axη * dηdx[iel,i,j,k] + axζ * dζdx[iel,i,j,k]) +
                        (ayξ * dξdy[iel,i,j,k] + ayη * dηdy[iel,i,j,k] + ayζ * dζdy[iel,i,j,k]) +
                        (azξ * dξdz[iel,i,j,k] + azη * dηdz[iel,i,j,k] + azζ * dζdz[iel,i,j,k]))
        end
    end
    return nothing
end

#-----------------------------------------------------------------------------
# Assembly tails -- the SAME DSS, halo exchange and mass division the five-field
# path uses, applied to a 3-wide and a 1-wide array instead of a 5-wide one.
#-----------------------------------------------------------------------------
@inline function _schur_assemble!(R, rel, nq::Int, params, op, vaux)
    mesh  = params.mesh
    npoin = Int(mesh.npoin)
    fill!(R, 0.0)
    DSS_rhs!(R, rel, mesh.connijk, Int(mesh.nelem), Int(mesh.ngl), nq,
             params.SD, params.AD)
    op.dss_cache !== nothing && assemble_mpi!(R, op.dss_cache)
    for q = 1:nq
        divide_by_mass_matrix!(@view(R[:, q]), vaux, params.Minv, nq, npoin, params.AD)
    end
    return R
end

"""
    schur_grad_fast!(gx, gy, gz, P, params, op, kern) -> nothing

`(gx, gy, gz) = Grad[P]` as the momentum rows of the operator see it, without
applying the operator.
"""
function schur_grad_fast!(gx, gy, gz, P, params, op, kern::SchurKernel)
    mesh = params.mesh
    met  = params.metrics
    _schur_grad_elements!(kern.rel3, P, kern.Pel, mesh.connijk,
                          params.basis.dψ, params.ω, met.Je,
                          met.dξdx, met.dξdy, met.dξdz,
                          met.dηdx, met.dηdy, met.dηdz,
                          met.dζdx, met.dζdy, met.dζdz,
                          Int(mesh.nelem), Int(mesh.ngl))
    R = _schur_assemble!(kern.R3, kern.rel3, 3, params, op, kern.vaux)
    @inbounds for ip in eachindex(P)
        gx[ip] = R[ip, 1]; gy[ip] = R[ip, 2]; gz[ip] = R[ip, 3]
    end
    return nothing
end

"""
    schur_divW_fast!(out, vx, vy, vz, params, op, kern) -> nothing

`out = Div[W v]`, the free-slip-masked divergence, without applying the
operator.
"""
function schur_divW_fast!(out, vx, vy, vz, params, op, kern::SchurKernel)
    mesh = params.mesh
    met  = params.metrics
    _schur_div_elements!(kern.rel1, vx, vy, vz, op.wall, op.wallx, op.wally,
                         kern.Mx, kern.My, kern.Mz, mesh.connijk,
                         params.basis.dψ, params.ω, met.Je,
                         met.dξdx, met.dξdy, met.dξdz,
                         met.dηdx, met.dηdy, met.dηdz,
                         met.dζdx, met.dζdy, met.dζdz,
                         Int(mesh.nelem), Int(mesh.ngl))
    R = _schur_assemble!(kern.R1, kern.rel1, 1, params, op, kern.vaux)
    @inbounds for ip in eachindex(out)
        out[ip] = R[ip, 1]
    end
    return nothing
end
