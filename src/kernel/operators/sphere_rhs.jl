#---------------------------------------------------------------------------------
# sphere_rhs.jl
#
# The right-hand side of the shallow water equations on the spherical shell,
# assembled by the SPECTRAL ELEMENT METHOD — no finite differences anywhere.
#
# Continuous Galerkin, strong form, exactly the convention of
# _expansion_inviscid!(…, NSD_2D, ContGal) in src/kernel/operators/rhs.jl, with
# the third flux component and the third metric direction added:
#
#   rhs_el[e,i,j,eq] -= ω_i ω_j J [ (∂F/∂x + ∂G/∂y + ∂H/∂z) - S ]
#
#   ∂F/∂x = dξdx ∂F/∂ξ + dηdx ∂F/∂η          (and likewise for G,y and H,z)
#   ∂F/∂ξ = Σ_k dψ[k,i] F[k,j]               the LGL differentiation matrix
#
# then DIRECT STIFFNESS SUMMATION (the `+=` into the global array, which is what
# makes the solution continuous across elements) and a division by the diagonal
# mass matrix:
#
#   ∂q/∂t = M⁻¹ Σ_e rhs_el .
#
# There is NO boundary term. The shell is a closed manifold: ∫_Γ vanishes
# identically, as the paper notes for the CG case (Marras, Kopera & Giraldo
# 2015, section 4.1, below Eq. 14). The "periodicity" is the topology itself,
# enforced by the node numbering built in sphere_mesh.jl.
#
# The fluxes come from the case's user_flux!, the sources from its user_source!,
# through the ordinary Jexpresso callbacks — this file never hard-codes an
# equation.
#
# ALSO HERE: the modal filter. The paper applies a Boyd-Vandeven filter at every
# step, "chosen to reduce by 5% the highest modes only", because the inviscid
# high-order solution of an advective problem grows spurious grid-scale modes.
# What is implemented is the exponential modal filter, which serves the same
# purpose and has the same knobs; see sphere_filter! below.
#
# S. Marras & contributors
#---------------------------------------------------------------------------------

export St_sphere_params
export build_sphere_params
export sphere_rhs!
export sphere_filter!
export build_sphere_filter_matrix
export sphere_relative_vorticity!


#---------------------------------------------------------------------------------
# Scratch space, allocated once and reused every stage of every step.
#---------------------------------------------------------------------------------
mutable struct St_sphere_params{TFloat}
    neqs ::Int
    F    ::Array{TFloat, 3}    # ngl × ngl × neqs
    G    ::Array{TFloat, 3}
    H    ::Array{TFloat, 3}
    S    ::Array{TFloat, 3}
    Filt ::Array{TFloat, 2}    # ngl × ngl modal filter (identity when disabled)
    qf   ::Array{TFloat, 3}    # ngl × ngl × neqs filter scratch
    acc  ::Array{TFloat, 2}    # npoin × neqs   filter accumulator
    lfilter ::Bool
end


function build_sphere_params(mesh::St_mesh, metrics::St_sphere_metrics, inputs;
                             neqs = 4, TF = TFloat)

    ngl = Int(mesh.ngl)

    lfilter = get(inputs, :lfilter, false) == true
    Filt    = lfilter ? build_sphere_filter_matrix(metrics.ξ, inputs, TF) :
                        Matrix{TF}(I, ngl, ngl)

    return St_sphere_params{TF}(neqs,
                                zeros(TF, ngl, ngl, neqs),
                                zeros(TF, ngl, ngl, neqs),
                                zeros(TF, ngl, ngl, neqs),
                                zeros(TF, ngl, ngl, neqs),
                                Filt,
                                zeros(TF, ngl, ngl, neqs),
                                zeros(TF, Int(mesh.npoin), neqs),
                                lfilter)
end


#---------------------------------------------------------------------------------
# sphere_rhs!(RHS, q, qe, mesh, metrics, sp, SVT)
#
# RHS  npoin × neqs, overwritten with ∂q/∂t.
# q    npoin × (neqs+1) conservative state (only the first neqs columns are read)
# qe   npoin × (neqs+1) reference state handed to the user callbacks
#---------------------------------------------------------------------------------
#
# TYPED FUNCTION BARRIER — the same pattern as _viscous_rhs_el_2d! in
# src/kernel/operators/rhs.jl, and it is not optional here.
#
# St_mesh is a Base.@kwdef struct whose fields carry NO type annotation, so
# `mesh.x`, `mesh.connijk`, `mesh.npoin` … are all `::Any`. Reading them inside
# the element loop makes every access a dynamic lookup that boxes its result:
# on the Galewsky run that cost 63 G allocations / 1.1 TiB and roughly 10x the
# runtime. Hoisting them into concretely-typed arguments once per call lets
# Julia specialise the kernel below, and the inner loops become allocation-free.
#
# Anything read per element or per node belongs in the signature. Anything read
# once per call (SD, the user callbacks' `mesh` argument) can stay dynamic.
#
function sphere_rhs!(RHS, q, qe,
                     mesh::St_mesh,
                     metrics::St_sphere_metrics,
                     sp::St_sphere_params,
                     SVT)

    _sphere_rhs_kernel!(RHS, q, qe,
                        mesh.connijk, mesh.x, mesh.y, mesh.z,
                        Int(mesh.nelem), Int(mesh.ngl), Int(mesh.npoin),
                        mesh, mesh.SD,
                        metrics.Je, metrics.Minv,
                        metrics.dξdx, metrics.dξdy, metrics.dξdz,
                        metrics.dηdx, metrics.dηdy, metrics.dηdz,
                        metrics.dψ, metrics.ω,
                        sp.F, sp.G, sp.H, sp.S, Int(sp.neqs), SVT)
    return RHS
end

function _sphere_rhs_kernel!(RHS::AbstractMatrix{TF}, q, qe,
                             connijk, xn::AbstractVector{TF},
                             yn::AbstractVector{TF}, zn::AbstractVector{TF},
                             nelem::Int, ngl::Int, npoin::Int,
                             mesh, SD,
                             Je, Minv::AbstractVector{TF},
                             dξdx, dξdy, dξdz,
                             dηdx, dηdy, dηdz,
                             dψ, ω,
                             F, G, H, S, neqs::Int, SVT) where {TF}

    fill!(RHS, zero(TF))

    @inbounds for iel = 1:nelem

        #--- fluxes and sources at the element's LGL nodes
        for j = 1:ngl, i = 1:ngl
            ip = connijk[iel,i,j]

            user_flux!(@view(F[i,j,:]), @view(G[i,j,:]), @view(H[i,j,:]),
                       SD, @view(q[ip,:]), @view(qe[ip,:]),
                       mesh, CL(), SVT; neqs = neqs, ip = ip)

            user_source!(@view(S[i,j,:]), @view(q[ip,:]), @view(qe[ip,:]),
                         npoin, CL(), SVT;
                         neqs = neqs,
                         x = xn[ip], y = yn[ip], z = zn[ip])
        end

        #--- surface divergence + direct stiffness summation
        for ieq = 1:neqs
            for j = 1:ngl
                ωj = ω[j]
                for i = 1:ngl

                    Jij  = Je[iel,i,j]
                    ωJac = ω[i]*ωj*Jij

                    dFdξ = dFdη = zero(TF)
                    dGdξ = dGdη = zero(TF)
                    dHdξ = dHdη = zero(TF)
                    for k = 1:ngl
                        dξ_ki = dψ[k,i]
                        dη_kj = dψ[k,j]
                        dFdξ += dξ_ki*F[k,j,ieq];  dFdη += dη_kj*F[i,k,ieq]
                        dGdξ += dξ_ki*G[k,j,ieq];  dGdη += dη_kj*G[i,k,ieq]
                        dHdξ += dξ_ki*H[k,j,ieq];  dHdη += dη_kj*H[i,k,ieq]
                    end

                    dFdx = dFdξ*dξdx[iel,i,j] + dFdη*dηdx[iel,i,j]
                    dGdy = dGdξ*dξdy[iel,i,j] + dGdη*dηdy[iel,i,j]
                    dHdz = dHdξ*dξdz[iel,i,j] + dHdη*dηdz[iel,i,j]

                    RHS[connijk[iel,i,j], ieq] -= ωJac*((dFdx + dGdy + dHdz) - S[i,j,ieq])
                end
            end
        end
    end

    #--- M⁻¹ : the mass matrix is diagonal under LGL (inexact) integration
    @inbounds for ieq = 1:neqs, ip = 1:npoin
        RHS[ip,ieq] *= Minv[ip]
    end

    return RHS
end


#---------------------------------------------------------------------------------
# build_sphere_filter_matrix(ξ, inputs, TF)
#
# The exponential modal filter, in nodal form: Filt = V Σ V⁻¹ where V is the
# Legendre Vandermonde on the LGL points and Σ damps the high modes,
#
#   σ_k = 1                                        for k ≤ k₀
#   σ_k = 1 - α ((k - k₀)/(N - k₀))^s              for k > k₀
#
# so σ_N = 1 - α: :filter_alpha is literally "how much the TOP mode is reduced
# by". The paper's Boyd-Vandeven filter "chosen to reduce by 5% the highest
# modes only" corresponds to α = 0.05 with a large s and k₀ near N, which is the
# default here.
#
#   :filter_alpha  damping of the highest mode          (default 0.05)
#   :filter_order  exponent s: larger = sharper cutoff  (default 8)
#   :filter_kcut   first damped mode k₀, as a fraction of N (default 2/3)
#---------------------------------------------------------------------------------
function build_sphere_filter_matrix(ξ, inputs, TF = TFloat)

    ngl = length(ξ)
    N   = ngl - 1

    α  = TF(get(inputs, :filter_alpha, 0.05))
    s  = TF(get(inputs, :filter_order, 8))
    k0 = max(0, min(N-1, round(Int, TF(get(inputs, :filter_kcut, 2/3))*N)))

    # Legendre Vandermonde: V[i,k] = P_{k-1}(ξ_i)
    V = zeros(TF, ngl, ngl)
    for i = 1:ngl
        p0 = one(TF); p1 = ξ[i]
        V[i,1] = p0
        ngl >= 2 && (V[i,2] = p1)
        for k = 2:N
            p0, p1 = p1, ((2k-1)*ξ[i]*p1 - (k-1)*p0)/k
            V[i,k+1] = p1
        end
    end

    σ = ones(TF, ngl)
    for k = k0+1:N
        σ[k+1] = one(TF) - α*((TF(k - k0)/TF(N - k0))^s)
    end

    return V * Diagonal(σ) * inv(V)
end


#---------------------------------------------------------------------------------
# sphere_filter!(q, mesh, metrics, sp)
#
# Apply the modal filter element by element (tensor product: rows then columns),
# then restore continuity with a MASS-WEIGHTED direct stiffness average,
#
#   q[ip] = Σ_e ω_i ω_j J q̃ / M[ip] ,
#
# which is the L² projection back onto the continuous space and therefore
# conserves ∫q exactly — the filter must not create or destroy mass.
#---------------------------------------------------------------------------------
function sphere_filter!(q, mesh::St_mesh, metrics::St_sphere_metrics, sp::St_sphere_params)

    sp.lfilter || return q

    # Typed function barrier; see the note above _sphere_rhs_kernel!.
    _sphere_filter_kernel!(q, mesh.connijk, Int(mesh.nelem), Int(mesh.ngl), Int(mesh.npoin),
                           metrics.Je, metrics.Minv, metrics.ω,
                           sp.Filt, sp.acc, sp.qf, sp.S, Int(sp.neqs))
    return q
end

function _sphere_filter_kernel!(q::AbstractMatrix{TF}, connijk,
                                nelem::Int, ngl::Int, npoin::Int,
                                Je, Minv::AbstractVector{TF}, ω,
                                Fm, acc, qf, S, neqs::Int) where {TF}

    fill!(acc, zero(TF))

    @inbounds for iel = 1:nelem

        # filter along ξ, then along η
        for ieq = 1:neqs
            for j = 1:ngl, i = 1:ngl
                t = zero(TF)
                for k = 1:ngl
                    t += Fm[i,k]*q[connijk[iel,k,j], ieq]
                end
                qf[i,j,ieq] = t
            end
        end
        for ieq = 1:neqs
            for i = 1:ngl
                for j = 1:ngl
                    t = zero(TF)
                    for k = 1:ngl
                        t += Fm[j,k]*qf[i,k,ieq]
                    end
                    S[i,j,ieq] = t             # reuse S as scratch
                end
            end
        end

        # mass-weighted accumulation
        for ieq = 1:neqs, j = 1:ngl, i = 1:ngl
            acc[connijk[iel,i,j], ieq] += ω[i]*ω[j]*Je[iel,i,j]*S[i,j,ieq]
        end
    end

    @inbounds for ieq = 1:neqs, ip = 1:npoin
        q[ip,ieq] = acc[ip,ieq]*Minv[ip]
    end

    return q
end


#---------------------------------------------------------------------------------
# sphere_relative_vorticity!(ζ, q, mesh, metrics)
#
# The RELATIVE VORTICITY of the tangential flow,
#
#   ζ = n̂ · (∇ₛ × u) ,   ∇ₛ × u = a¹ × ∂u/∂ξ + a² × ∂u/∂η ,
#
# with u = φu/φ the Cartesian velocity. This is the field the Galewsky test is
# judged on — the paper (and Galewsky et al.) plot relative vorticity at day 6,
# where the barotropic instability has rolled the jet into a train of vortices.
# The height field barely moves by comparison, so h alone makes a developing
# instability almost invisible.
#
# Computed element-wise from the metric terms and then made continuous with the
# same mass-weighted direct stiffness average used by the filter, so ζ is a
# proper nodal field on the shell rather than a discontinuous per-element one.
#---------------------------------------------------------------------------------
function sphere_relative_vorticity!(ζ, q, mesh::St_mesh,
                                    metrics::St_sphere_metrics, sp::St_sphere_params)

    ngl = Int(mesh.ngl)
    dψ  = metrics.dψ
    ω   = metrics.ω

    acc = sp.acc                      # npoin × neqs scratch (column 1 only)
    fill!(acc, zero(eltype(acc)))

    @inbounds for iel = 1:mesh.nelem
        for j = 1:ngl, i = 1:ngl

            # ∂u/∂ξ and ∂u/∂η of the CARTESIAN velocity
            dξx = dξy = dξz = 0.0
            dηx = dηy = dηz = 0.0
            for k = 1:ngl
                ipk = mesh.connijk[iel,k,j]
                φk  = q[ipk,1]
                dξx += dψ[k,i]*q[ipk,2]/φk; dξy += dψ[k,i]*q[ipk,3]/φk; dξz += dψ[k,i]*q[ipk,4]/φk

                ipl = mesh.connijk[iel,i,k]
                φl  = q[ipl,1]
                dηx += dψ[k,j]*q[ipl,2]/φl; dηy += dψ[k,j]*q[ipl,3]/φl; dηz += dψ[k,j]*q[ipl,4]/φl
            end

            a1x, a1y, a1z = metrics.dξdx[iel,i,j], metrics.dξdy[iel,i,j], metrics.dξdz[iel,i,j]
            a2x, a2y, a2z = metrics.dηdx[iel,i,j], metrics.dηdy[iel,i,j], metrics.dηdz[iel,i,j]

            cx = (a1y*dξz - a1z*dξy) + (a2y*dηz - a2z*dηy)
            cy = (a1z*dξx - a1x*dξz) + (a2z*dηx - a2x*dηz)
            cz = (a1x*dξy - a1y*dξx) + (a2x*dηy - a2y*dηx)

            zt = cx*metrics.nx[iel,i,j] + cy*metrics.ny[iel,i,j] + cz*metrics.nz[iel,i,j]

            acc[mesh.connijk[iel,i,j], 1] += ω[i]*ω[j]*metrics.Je[iel,i,j]*zt
        end
    end

    @inbounds for ip = 1:mesh.npoin
        ζ[ip] = acc[ip,1]*metrics.Minv[ip]
    end

    return ζ
end
