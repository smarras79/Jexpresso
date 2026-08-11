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
function sphere_rhs!(RHS, q, qe,
                     mesh::St_mesh,
                     metrics::St_sphere_metrics,
                     sp::St_sphere_params,
                     SVT)

    ngl  = Int(mesh.ngl)
    neqs = sp.neqs
    dψ   = metrics.dψ
    ω    = metrics.ω

    F, G, H, S = sp.F, sp.G, sp.H, sp.S

    fill!(RHS, zero(eltype(RHS)))

    @inbounds for iel = 1:mesh.nelem

        #--- fluxes and sources at the element's LGL nodes
        for j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[iel,i,j]

            user_flux!(@view(F[i,j,:]), @view(G[i,j,:]), @view(H[i,j,:]),
                       mesh.SD, @view(q[ip,:]), @view(qe[ip,:]),
                       mesh, CL(), SVT; neqs = neqs, ip = ip)

            user_source!(@view(S[i,j,:]), @view(q[ip,:]), @view(qe[ip,:]),
                         mesh.npoin, CL(), SVT;
                         neqs = neqs,
                         x = mesh.x[ip], y = mesh.y[ip], z = mesh.z[ip])
        end

        #--- surface divergence + direct stiffness summation
        for ieq = 1:neqs
            for j = 1:ngl
                ωj = ω[j]
                for i = 1:ngl

                    Jij  = metrics.Je[iel,i,j]
                    ωJac = ω[i]*ωj*Jij

                    dFdξ = dFdη = 0.0
                    dGdξ = dGdη = 0.0
                    dHdξ = dHdη = 0.0
                    for k = 1:ngl
                        dξ_ki = dψ[k,i]
                        dη_kj = dψ[k,j]
                        dFdξ += dξ_ki*F[k,j,ieq];  dFdη += dη_kj*F[i,k,ieq]
                        dGdξ += dξ_ki*G[k,j,ieq];  dGdη += dη_kj*G[i,k,ieq]
                        dHdξ += dξ_ki*H[k,j,ieq];  dHdη += dη_kj*H[i,k,ieq]
                    end

                    dFdx = dFdξ*metrics.dξdx[iel,i,j] + dFdη*metrics.dηdx[iel,i,j]
                    dGdy = dGdξ*metrics.dξdy[iel,i,j] + dGdη*metrics.dηdy[iel,i,j]
                    dHdz = dHdξ*metrics.dξdz[iel,i,j] + dHdη*metrics.dηdz[iel,i,j]

                    RHS[mesh.connijk[iel,i,j], ieq] -= ωJac*((dFdx + dGdy + dHdz) - S[i,j,ieq])
                end
            end
        end
    end

    #--- M⁻¹ : the mass matrix is diagonal under LGL (inexact) integration
    @inbounds for ieq = 1:neqs, ip = 1:mesh.npoin
        RHS[ip,ieq] *= metrics.Minv[ip]
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

    ngl  = Int(mesh.ngl)
    neqs = sp.neqs
    Fm   = sp.Filt
    ω    = metrics.ω
    acc  = sp.acc
    qf   = sp.qf

    fill!(acc, zero(eltype(acc)))

    @inbounds for iel = 1:mesh.nelem

        # filter along ξ, then along η
        for ieq = 1:neqs
            for j = 1:ngl, i = 1:ngl
                t = 0.0
                for k = 1:ngl
                    t += Fm[i,k]*q[mesh.connijk[iel,k,j], ieq]
                end
                qf[i,j,ieq] = t
            end
        end
        for ieq = 1:neqs
            for i = 1:ngl
                for j = 1:ngl
                    t = 0.0
                    for k = 1:ngl
                        t += Fm[j,k]*qf[i,k,ieq]
                    end
                    sp.S[i,j,ieq] = t          # reuse S as scratch
                end
            end
        end

        # mass-weighted accumulation
        for ieq = 1:neqs, j = 1:ngl, i = 1:ngl
            acc[mesh.connijk[iel,i,j], ieq] += ω[i]*ω[j]*metrics.Je[iel,i,j]*sp.S[i,j,ieq]
        end
    end

    @inbounds for ieq = 1:neqs, ip = 1:mesh.npoin
        q[ip,ieq] = acc[ip,ieq]*metrics.Minv[ip]
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
