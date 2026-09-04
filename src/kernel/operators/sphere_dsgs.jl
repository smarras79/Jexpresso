#---------------------------------------------------------------------------------
# sphere_dsgs.jl
#
# DYNAMIC SUB-GRID-SCALE (residual-based) ARTIFICIAL VISCOSITY on the spherical
# shell — the DynSGS of
#
#   Marras, S., Nazarov, M. & Giraldo, F. X. (2015), "Stabilized high-order
#   Galerkin methods based on a parameter-free dynamic SGS model for LES",
#   J. Comput. Phys. 301, 77-101,
#
# in the form the flat cases run it (:visc_model => DSGS(), see
# src/kernel/physics/SGS.jl), carried to the shell's own RHS.
#
# WHAT IT REPLACES. The shell stabilises with a CONSTANT ν∇ₛ²(φu) and the
# modal filter (sphere_rhs.jl). Both act everywhere at all times, whether or
# not anything needs damping: the forced-turbulence deck of
# SWsphere_ScottPolvani has to keep ν small so as not to eat the forced band,
# and then leans on the filter for the grid scale. DynSGS puts viscosity only
# where the discrete equations are NOT being satisfied — where the residual is
# large, i.e. at under-resolved scales — and none elsewhere.
#
# THE COEFFICIENT, per element e (Marras et al. Eqs. 8-10):
#
#   ν_res|e = C₁ Δ² max_i ‖R_i‖∞,e / ‖q_i − ⟨q_i⟩‖∞,Ω
#   ν_max|e = C₂ Δ  ‖ |u| + √φ ‖∞,e
#   ν|e     = max(0, min(ν_max|e, ν_res|e))
#
#   Δ       = (shortest corner-to-corner edge of the element)/ngl, the LGL
#             node spacing the element actually resolves;
#   R_i     = ∂q_i/∂t − RHS_i(q), the strong-form residual of equation i, with
#             ∂q_i/∂t the BDF2 stencil (3qⁿ − 4qⁿ⁻¹ + qⁿ⁻²)/(2Δt) over the last
#             two COMPLETED steps and RHS the assembled, M⁻¹-scaled inviscid
#             right-hand side of the current stage (forcing and Rayleigh
#             friction included — everything but the viscous term itself);
#   ⟨q_i⟩, ‖q_i − ⟨q_i⟩‖∞,Ω  the area mean and the ∞-norm of the deviation
#             from it over the domain — rank-local unless :ldsgs_global_norms,
#             see _dsgs_norm_scope in SGS.jl for why local is the default;
#   C₁ = 1, C₂ = 0.5 the paper's constants (:dsgs_C1, :dsgs_C2).
#
# The cap ν_max is what makes the model explicit-stable by construction:
# Δt ν/Δ² ≤ C₂ (|u|+c) Δt/Δ ≈ C₂·CFL.
#
# FLOORS ON THE NORMALISATION. The paper divides by ‖q − ⟨q⟩‖∞,Ω. A resting
# layer of uniform depth — the Scott & Polvani start — has that deviation
# EXACTLY zero for every equation, and the residual over zero would pin ν at
# the cap everywhere until the flow has developed, which is maximal
# dissipation precisely when nothing needs it. As on the flat θ path (the
# momentum floor 1e-3·ρc in SGS.jl) each denominator gets a floor at a small
# fraction of its natural scale, φ̄ for the geopotential and φ̄√φ̄ for the
# momentum (:dsgs_floor, default 1e-3). Once the actual deviations exceed the
# floor it plays no role.
#
# WHAT IS DIFFUSED. The shell's viscous operator is ν∇ₛ²q on the CONSERVED
# Cartesian momentum (φu, φv, φw), the δν∇²(φu) of Marras, Kopera & Giraldo
# (2015) Eq. (8b), so ν is a KINEMATIC coefficient [m²/s] and the element
# value above is applied as it is, times the per-equation multiplier :μ
# ([0, 1, 1, 1]: no mass diffusion, Marras Eq. 10). The Laplace-Beltrami weak
# form and the direct stiffness summation are the ones of sphere_rhs.jl; the
# only change there is that the coefficient is ν[ieq]·ν_e[iel] instead of
# ν[ieq].
#
# TWO ASSEMBLY PASSES. The residual needs the inviscid RHS BEFORE the viscous
# term is added, so on a DynSGS run sphere_rhs! assembles and M⁻¹-scales the
# inviscid part first, computes ν_e from it, and only then assembles the
# viscous part into scratch and adds it. The constant-ν path keeps its single
# pass and is bit-for-bit what it was.
#
# HISTORY. qⁿ and qⁿ⁻¹ are the states after the last two completed steps,
# rolled by the step limiter of sphere_time_loop.jl exactly once per step
# (the per-stage states of an RK scheme are no use to a time derivative).
# Both start at q⁰, so the first residual is zero and the run starts without
# viscosity, as the flat implementation does.
#
# DECK SWITCHES
#   :lvisc               => true
#   :visc_model          => DSGS()
#   :μ                   => [0.0, 1.0, 1.0, 1.0]   per-equation multiplier
#   :dsgs_C1, :dsgs_C2   => 1.0, 0.5
#   :dsgs_floor          => 1.0e-3
#   :ldsgs_global_norms  => false
#
# S. Marras & contributors
#---------------------------------------------------------------------------------

export St_sphere_dsgs
export build_sphere_dsgs
export sphere_dsgs_viscosity!
export sphere_dsgs_roll!
export sphere_dsgs_nodal!
export sphere_dsgs_diagnostics
export sphere_dsgs_cap_estimate


mutable struct St_sphere_dsgs{TF}
    C1      ::TF
    C2      ::TF
    floor   ::TF
    lglobal ::Bool
    Δ       ::Vector{TF}     # nelem : Marras's Δ = (shortest corner edge)/ngl
    νel     ::Vector{TF}     # nelem : ν_e of the last evaluation [m²/s]
    νnode   ::Vector{TF}     # npoin : max over the elements at a node (output only)
    q1      ::Matrix{TF}     # npoin × ncol : qⁿ    (after the last completed step)
    q2      ::Matrix{TF}     # npoin × ncol : qⁿ⁻¹
    #--- diagnostics of the last evaluation (this rank)
    νmax    ::TF
    νsum    ::TF             # Σ_e ν_e, for a mean
    ncap    ::Int            # elements sitting at the cap ν_max
end


#
# Per-element Δ: the shortest of the four corner-to-corner edges, divided by
# ngl — Marras's min(Δx, Δy)/(N+1). From the 3-D coordinates, since
# mesh.Δelem is computed on (x, y) by the flat path and is meaningless on a
# shell.
#
function _sphere_element_delta(mesh::St_mesh, TF)
    ngl   = Int(mesh.ngl)
    nelem = Int(mesh.nelem)
    crd   = mesh.coords
    Δ = zeros(TF, nelem)
    @inbounds for iel = 1:nelem
        c1 = mesh.connijk[iel, 1,   1  ]
        c2 = mesh.connijk[iel, ngl, 1  ]
        c3 = mesh.connijk[iel, ngl, ngl]
        c4 = mesh.connijk[iel, 1,   ngl]
        d(a, b) = sqrt((crd[1,a]-crd[1,b])^2 + (crd[2,a]-crd[2,b])^2 + (crd[3,a]-crd[3,b])^2)
        Δ[iel] = min(d(c1,c2), d(c2,c3), d(c3,c4), d(c4,c1)) / ngl
    end
    return Δ
end


function build_sphere_dsgs(mesh::St_mesh, inputs, q0::AbstractMatrix; neqs = 4,
                           TF = TFloat, verbose = true)

    get(inputs, :lvisc, false) == true || return nothing
    get(inputs, :visc_model, nothing) == DSGS() || return nothing

    neqs >= 4 || error(" # ERROR sphere_dsgs.jl: the shell DynSGS needs the 4-equation state [φ, φu, φv, φw].")

    C1    = TF(get(inputs, :dsgs_C1, 1.0))
    C2    = TF(get(inputs, :dsgs_C2, 0.5))
    floor = TF(get(inputs, :dsgs_floor, 1.0e-3))
    lglob = get(inputs, :ldsgs_global_norms, false) == true
    C1 >= 0 && C2 >= 0 && floor >= 0 ||
        error(" # ERROR sphere_dsgs.jl: :dsgs_C1, :dsgs_C2 and :dsgs_floor must be ≥ 0.")

    npoin = Int(mesh.npoin)
    Δ     = _sphere_element_delta(mesh, TF)

    ds = St_sphere_dsgs{TF}(C1, C2, floor, lglob, Δ,
                            zeros(TF, Int(mesh.nelem)), zeros(TF, npoin),
                            Matrix{TF}(q0[1:npoin, :]), Matrix{TF}(q0[1:npoin, :]),
                            zero(TF), zero(TF), 0)

    # Δ range over the WHOLE grid for the printout. `init`: a rank can own no
    # element at all (the x-y box partition of a sphere leaves the corner boxes
    # of the bounding square empty at large rank counts), and a bare
    # maximum() over an empty vector throws.
    comm = get_mpi_comm()
    Δlo  = minimum(Δ; init = TF(Inf)); Δhi = maximum(Δ; init = zero(TF))
    if MPI.Comm_size(comm) > 1
        Δlo = MPI.Allreduce(Δlo, MPI.MIN, comm); Δhi = MPI.Allreduce(Δhi, MPI.MAX, comm)
    end
    if verbose
        println(" # SHELL DynSGS (Marras, Nazarov & Giraldo 2015) ..................")
        @printf(" #   ν_e = max(0, min(C₂ Δ (|u|+√φ), C₁ Δ² max_i ‖R_i‖/‖q_i-⟨q_i⟩‖)),  C₁ = %.3g, C₂ = %.3g\n", C1, C2)
        @printf(" #   Δ = shortest corner edge / ngl ∈ [%.4e, %.4e] m ; normalisation floor %.1e of the natural scales ; norms %s\n",
                Δlo, Δhi, floor, lglob ? "global" : "rank-local")
        println(" # SHELL DynSGS ................................................. DONE")
    end
    return ds
end


#
# Roll the BDF2 history after a completed step: qⁿ⁻¹ ← qⁿ, qⁿ ← u.
#
function sphere_dsgs_roll!(ds::St_sphere_dsgs{TF}, u) where {TF}
    npoin = size(ds.q1, 1)
    ncol  = size(ds.q1, 2)
    @inbounds for j = 1:ncol, ip = 1:npoin
        ds.q2[ip, j] = ds.q1[ip, j]
        ds.q1[ip, j] = u[ip, j]
    end
    return ds
end


#---------------------------------------------------------------------------------
# sphere_dsgs_viscosity!(ds, q, RHS, mesh, metrics, Δt, neqs)
#
# ν_e for every element of this rank, from the current state q and the
# assembled, M⁻¹-scaled inviscid RHS (npoin × neqs). Writes ds.νel.
#---------------------------------------------------------------------------------
function sphere_dsgs_viscosity!(ds::St_sphere_dsgs{TF}, q, RHS::AbstractMatrix{TF},
                                mesh::St_mesh, metrics::St_sphere_metrics,
                                Δt, neqs::Int) where {TF}
    _sphere_dsgs_kernel!(ds, q, RHS, mesh.connijk, metrics.M, mesh.gip2owner,
                         MPI.Comm_rank(get_mpi_comm()),
                         Int(mesh.nelem), Int(mesh.ngl), Int(mesh.npoin), TF(Δt))
    return ds
end

function _sphere_dsgs_kernel!(ds::St_sphere_dsgs{TF}, q::AbstractMatrix{TF},
                              RHS::AbstractMatrix{TF}, connijk,
                              M::AbstractVector{TF}, gip2owner, rank::Int,
                              nelem::Int, ngl::Int, npoin::Int, Δt::TF) where {TF}

    comm = get_mpi_comm()
    eps  = TF(1.0e-16)

    #--- pass 1: area means ⟨q_i⟩ over this rank's owned nodes (or the domain)
    s1 = s2 = s3 = s4 = sm = zero(TF)
    @inbounds for ip = 1:npoin
        gip2owner[ip] == rank || continue
        w   = M[ip]
        s1 += w*q[ip,1]; s2 += w*q[ip,2]; s3 += w*q[ip,3]; s4 += w*q[ip,4]
        sm += w
    end
    if ds.lglobal && MPI.Comm_size(comm) > 1
        sums = TF[s1, s2, s3, s4, sm]
        MPI.Allreduce!(sums, MPI.SUM, comm)
        s1, s2, s3, s4, sm = sums[1], sums[2], sums[3], sums[4], sums[5]
    end
    a1 = s1/sm; a2 = s2/sm; a3 = s3/sm; a4 = s4/sm

    #--- pass 2: ‖q_i − ⟨q_i⟩‖∞
    d1 = d2 = d3 = d4 = zero(TF)
    @inbounds for ip = 1:npoin
        d1 = max(d1, abs(q[ip,1] - a1)); d2 = max(d2, abs(q[ip,2] - a2))
        d3 = max(d3, abs(q[ip,3] - a3)); d4 = max(d4, abs(q[ip,4] - a4))
    end
    if ds.lglobal && MPI.Comm_size(comm) > 1
        mx = TF[d1, d2, d3, d4]
        MPI.Allreduce!(mx, MPI.MAX, comm)
        d1, d2, d3, d4 = mx[1], mx[2], mx[3], mx[4]
    end
    # floors: a fraction of the natural scales φ̄ and φ̄√φ̄ (see the header)
    φ̄  = abs(a1)
    fφ = ds.floor*φ̄
    fm = ds.floor*φ̄*sqrt(φ̄)
    d1 = max(d1, fφ, eps); d2 = max(d2, fm, eps); d3 = max(d3, fm, eps); d4 = max(d4, fm, eps)

    #--- pass 3: the element residual, the cap, ν_e
    C1 = ds.C1; C2 = ds.C2
    i2Δt = one(TF)/(2Δt)
    νmax = zero(TF); νsum = zero(TF); ncap = 0
    @inbounds for iel = 1:nelem
        n1 = n2 = n3 = n4 = zero(TF)
        umax = zero(TF)
        for j = 1:ngl, i = 1:ngl
            ip = connijk[iel, i, j]
            r1 = abs((3q[ip,1] - 4ds.q1[ip,1] + ds.q2[ip,1])*i2Δt - RHS[ip,1])
            r2 = abs((3q[ip,2] - 4ds.q1[ip,2] + ds.q2[ip,2])*i2Δt - RHS[ip,2])
            r3 = abs((3q[ip,3] - 4ds.q1[ip,3] + ds.q2[ip,3])*i2Δt - RHS[ip,3])
            r4 = abs((3q[ip,4] - 4ds.q1[ip,4] + ds.q2[ip,4])*i2Δt - RHS[ip,4])
            n1 = max(n1, r1); n2 = max(n2, r2); n3 = max(n3, r3); n4 = max(n4, r4)
            φ  = q[ip,1]
            umag = sqrt(q[ip,2]^2 + q[ip,3]^2 + q[ip,4]^2)/φ
            umax = max(umax, umag + sqrt(max(φ, zero(TF))))
        end
        Δ    = ds.Δ[iel]
        νres = C1*Δ*Δ*max(n1/d1, n2/d2, n3/d3, n4/d4)
        νcap = C2*Δ*umax
        ν    = max(zero(TF), min(νcap, νres))
        ds.νel[iel] = ν
        νmax  = max(νmax, ν); νsum += ν
        νres >= νcap && (ncap += 1)
    end
    ds.νmax = νmax; ds.νsum = νsum; ds.ncap = ncap
    return ds
end


#
# Nodal field for output: the largest ν_e among the elements sharing a node.
#
function sphere_dsgs_nodal!(ds::St_sphere_dsgs{TF}, mesh::St_mesh) where {TF}
    fill!(ds.νnode, zero(TF))
    ngl = Int(mesh.ngl)
    connijk = mesh.connijk
    @inbounds for iel = 1:Int(mesh.nelem)
        ν = ds.νel[iel]
        for j = 1:ngl, i = 1:ngl
            ip = connijk[iel, i, j]
            ds.νnode[ip] = max(ds.νnode[ip], ν)
        end
    end
    return ds.νnode
end


#
# Reduced diagnostics of the last evaluation: max ν_e, mean ν_e, fraction of
# elements at the cap. Collective.
#
function sphere_dsgs_diagnostics(ds::St_sphere_dsgs{TF}, mesh::St_mesh) where {TF}
    comm  = get_mpi_comm()
    nelem = Int(mesh.nelem)
    νmax  = ds.νmax; νsum = ds.νsum; ncap = ds.ncap; ne = nelem
    if MPI.Comm_size(comm) > 1
        νmax = MPI.Allreduce(νmax, MPI.MAX, comm)
        νsum = MPI.Allreduce(νsum, MPI.SUM, comm)
        ncap = MPI.Allreduce(ncap, MPI.SUM, comm)
        ne   = MPI.Allreduce(ne,   MPI.SUM, comm)
    end
    return (νmax = νmax, νmean = νsum/max(ne, 1), fcap = ncap/max(ne, 1))
end


#
# The largest ν the cap can ever hand out for the state q — C₂ max(Δ)
# max(|u|+√φ) — for the diffusive branch of the CFL condition. Collective.
#
function sphere_dsgs_cap_estimate(ds::St_sphere_dsgs{TF}, q, npoin::Int) where {TF}
    cmax = zero(TF)
    @inbounds for ip = 1:npoin
        φ = q[ip,1]
        cmax = max(cmax, sqrt(q[ip,2]^2 + q[ip,3]^2 + q[ip,4]^2)/φ + sqrt(max(φ, zero(TF))))
    end
    comm = get_mpi_comm()
    Δmax = maximum(ds.Δ; init = zero(TF))     # init: this rank may own no element
    if MPI.Comm_size(comm) > 1
        cmax = MPI.Allreduce(cmax, MPI.MAX, comm)
        Δmax = MPI.Allreduce(Δmax, MPI.MAX, comm)
    end
    return ds.C2*Δmax*cmax
end
