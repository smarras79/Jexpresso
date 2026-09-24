#-------------------------------------------------------------------------------
# Numerical schlieren from the density field.
#
# A schlieren photograph records the deflection of light by refractive-index
# gradients, and for a gas the refractive index is an affine function of
# density (Gladstone-Dale), so what the picture actually shows is |∇ρ|. The
# standard numerical analogue (Quirk, ICASE 94-75; Hadjadj & Kudryavtsev,
# Shock Waves 14 (2005) 169) maps that magnitude through a decaying
# exponential
#
#     S = exp( -k · |∇ρ| / max_Ω |∇ρ| )
#
# so that S ∈ [e^-k, 1]: 1 in uniform flow, dark where the density gradient
# is steep. The exponential is what makes weak features (an expansion fan, a
# slip line, the contact behind a Mach stem) visible next to a strong bow
# shock — a linear |∇ρ| map is saturated by the strongest shock and shows
# nothing else. k controls that contrast; 20 is the usual choice, and it is
# `inputs[:schlieren_k]` here.
#
# BOTH fields are written out, because they answer different questions:
#
#   schlieren_grad_rho  |∇ρ| in kg/m⁴ — quantitative, use it to locate or
#                       threshold a shock
#   schlieren           S, dimensionless in [e^-k, 1] — the picture; plot
#                       this one in greyscale, reversed, for the familiar
#                       schlieren look
#
# Enable per case with `:lschlieren => true` (see
# problems/CompEuler/ffs_step/user_inputs.jl). Off by default, and computed
# only at output times, so it costs nothing in the time loop.
#-------------------------------------------------------------------------------

const _SCHLIEREN_DIM_WARNED = Ref(false)

#
# ∇ρ on the SEM grid, then the two fields above, into schlieren[1:npoin, 1:2].
#
# The gradient is the collocation derivative inside each element,
#
#     ∂ρ/∂ξ|_(i,j) = Σ_m dψ[m,i]·ρ(m,j),   ∂ρ/∂η|_(i,j) = Σ_m dψ[m,j]·ρ(i,m)
#
# pushed to physical space with the element metrics — the same construction
# the viscous operator in rhs.jl uses, so the gradient seen here is the one
# the discretisation actually has.
#
# A CG node shared by several elements gets one value per element, and they
# differ at exactly the places this diagnostic is meant to highlight (a shock
# inside an element face). They are averaged by multiplicity rather than
# L²-projected: this is a visualisation field, and the cheap average removes
# the element-boundary striping that picking one contributing element leaves
# behind. CAVEAT for MPI runs: only rank-local contributions are averaged, so
# a node on a partition boundary is averaged over the elements this rank
# owns. The seam is invisible at plotting contrast but it is there.
#
function compute_schlieren!(schlieren::AbstractMatrix{TT},
                            uaux::AbstractMatrix{TT},
                            mesh, basis, metrics, ::NSD_2D;
                            k::Real = 20.0,
                            irho::Int = 1) where {TT<:AbstractFloat}

    npoin = Int(mesh.npoin)
    nelem = Int(mesh.nelem)
    ngl   = Int(mesh.ngl)
    connijk = mesh.connijk
    dψ      = basis.dψ

    kT = TT(k)

    fill!(schlieren, zero(TT))
    nshare = zeros(Int32, npoin)

    @inbounds for iel = 1:nelem
        for j = 1:ngl
            for i = 1:ngl
                ip = connijk[iel,i,j,1]

                dρdξ = zero(TT)
                dρdη = zero(TT)
                for m = 1:ngl
                    dρdξ += dψ[m,i]*uaux[connijk[iel,m,j,1], irho]
                    dρdη += dψ[m,j]*uaux[connijk[iel,i,m,1], irho]
                end

                dρdx = dρdξ*metrics.dξdx[iel,i,j] + dρdη*metrics.dηdx[iel,i,j]
                dρdy = dρdξ*metrics.dξdy[iel,i,j] + dρdη*metrics.dηdy[iel,i,j]

                schlieren[ip,1] += sqrt(dρdx*dρdx + dρdy*dρdy)
                nshare[ip]      += one(Int32)
            end
        end
    end

    @inbounds for ip = 1:npoin
        if nshare[ip] > 1
            schlieren[ip,1] /= nshare[ip]
        end
    end

    # max_Ω |∇ρ| and the density range are DOMAIN quantities, so they have to
    # be reduced across ranks — otherwise each partition would normalise
    # against its own maximum and the image would change brightness from one
    # partition to the next. One Allreduce carries all three (min via -max).
    red = TT[npoin > 0 ? maximum(@view schlieren[1:npoin, 1]) : zero(TT),
             npoin > 0 ? maximum(@view uaux[1:npoin, irho])   : typemin(TT),
             npoin > 0 ? -minimum(@view uaux[1:npoin, irho])  : typemin(TT)]
    MPI.Allreduce!(red, MPI.MAX, get_mpi_comm())
    gmax = red[1]
    ρmax = red[2]
    ρmin = -red[3]

    # Is there any density variation to see?
    #
    # Testing `gmax > 0` is NOT enough, and getting this wrong is loud: the
    # collocation derivative of an exactly uniform field is round-off, not
    # zero (≈ 4e-13 kg/m⁴ on this grid, from the metric products and the
    # ngl-term sum). Normalising by THAT stretches pure round-off across the
    # whole [e^-k, 1] range and paints a full-contrast speckle image of
    # nothing — which is exactly what the t = 0 frame of every
    # supersonic-inflow case would show.
    #
    # So the test is on the relative density RANGE, which is scale-free and
    # needs no length scale: below 1e-10 of the density itself there is no
    # signal, only round-off, and the honest picture is a blank one.
    ρscale = max(abs(ρmax), abs(ρmin))
    uniform = (gmax <= zero(TT)) ||
              (ρmax - ρmin <= TT(1.0e-10)*max(ρscale, eps(TT)))

    if uniform
        @inbounds for ip = 1:npoin
            schlieren[ip,1] = zero(TT)   # report the round-off as the 0 it is
            schlieren[ip,2] = one(TT)    # S = exp(0) = 1, a blank frame
        end
    else
        inv_gmax = one(TT)/gmax
        @inbounds for ip = 1:npoin
            schlieren[ip,2] = exp(-kT*schlieren[ip,1]*inv_gmax)
        end
    end

    return schlieren
end

#
# Fallback for 1D/3D: leave the buffer alone rather than write silent zeros
# that look like "no gradients anywhere". Warned once per session.
#
function compute_schlieren!(schlieren::AbstractMatrix, uaux::AbstractMatrix,
                            mesh, basis, metrics, SD; kwargs...)
    if !_SCHLIEREN_DIM_WARNED[]
        _SCHLIEREN_DIM_WARNED[] = true
        @warn "compute_schlieren! is implemented for NSD_2D only; :lschlieren ignored for $(typeof(SD))."
    end
    return nothing
end

#
# Output-time entry point. Returns the npoin×2 buffer, or `nothing` when the
# case has not asked for a schlieren — `nothing` is what write_vtk checks, so
# a case that leaves :lschlieren off pays nothing and gets no extra fields.
#
function maybe_compute_schlieren(inputs, p, u)

    get(inputs, :lschlieren, false) == true || return nothing
    p.SD isa NSD_2D                         || return compute_schlieren!(
                                                   zeros(eltype(p.uaux), 1, 2),
                                                   p.uaux, p.mesh, p.basis,
                                                   p.metrics, p.SD)

    npoin = Int(p.mesh.npoin)
    TT    = eltype(p.uaux)

    # Refresh uaux from the flat integrator state, exactly as write_vtk does
    # a moment later, so the density differentiated here is the density
    # written to the file and not the previous RHS evaluation's.
    u2uaux!(p.uaux, u, p.qp.neqs, npoin)

    schlieren = zeros(TT, npoin, 2)
    compute_schlieren!(schlieren, p.uaux, p.mesh, p.basis, p.metrics, p.SD;
                       k = get(inputs, :schlieren_k, 20.0))
    return schlieren
end
