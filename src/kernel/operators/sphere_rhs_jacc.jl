#---------------------------------------------------------------------------------
#
# sphere_rhs_jacc.jl — the SEM right-hand side on the CLOSED SPHERICAL SHELL,
#                      written once with JACC.jl and run unchanged on
#
#                          threads (multi-core CPU) | CUDA | AMDGPU | oneAPI | Metal
#
# The shell's counterpart to src/kernel/operators/rhs_jacc.jl, and it needs to be
# a file of its own because the shell does not go through `rhs!` AT ALL. It has
# its own right-hand side (sphere_rhs!, sphere_rhs.jl) and its own time loop
# (sphere_time_loop.jl), so the flat path's `:ljacc` switch never reaches it.
#
# WHAT IS DIFFERENT FROM THE FLAT 2D KERNELS
# ------------------------------------------
# One thing of substance: the metric is 3x2, not 2x2. The shell is a 2D manifold
# embedded in 3D, so the two contravariant basis vectors
#
#     a¹ = (dξdx, dξdy, dξdz)      a² = (dηdx, dηdy, dηdz)
#
# each carry a z component, the flux has a THIRD component H, and the surface
# divergence is
#
#     ∇ₛ·F = dFdx + dGdy + dHdz
#
# with each term assembled from both ξ and η. The viscous term is likewise the
# LAPLACE-BELTRAMI operator, not the flat Laplacian — see _sphere_visc_el! in
# sphere_rhs.jl, which this mirrors exactly.
#
# Three things are SIMPLER here than on the flat path:
#
#   * a closed shell has no boundary, so there is no boundary kernel at all and
#     no user_bc_dirichlet_gpu in the picture;
#   * the shell's ODE state is a MATRIX (npoin x neqs+1), not a flat vector, so
#     there is no u -> uaux repacking;
#   * the direct stiffness summation reuses _sphere_dss_scale! verbatim on the
#     host, which is what makes the parallel semantics identical to the serial
#     path by construction rather than by argument.
#
# THE FLUX TUPLE IS EQUATION-MAJOR
# --------------------------------
# The shell's user_flux_gpu returns 3*neqs values ordered
#
#     (F₁,G₁,H₁, F₂,G₂,H₂, …)
#
# i.e. the three components of equation ieq's flux vector are adjacent, at
# 3(ieq-1)+1 … 3(ieq-1)+3. That is NOT the flat path's convention, where the
# whole of F comes before the whole of G. It follows the host user_flux!(F, G, H,
# …), which writes F[ieq], G[ieq], H[ieq], and it is what the case file already
# returns — so it is the convention, and the kernels below unpack it.
#
# NO ATOMICS — same reason as the flat path. The element kernels write only their
# own slot of rhs_el and a node-parallel pass gathers through a fixed CSR map, so
# the answer is identical on every backend and at every thread count. See the
# header of rhs_jacc.jl for the full argument.
#
# TURNING IT ON
#     user_inputs.jl:  :ljacc => true          (with :backend left at CPU())
#     shell:           julia --project=. -t 8 …
#
#---------------------------------------------------------------------------------


#---------------------------------------------------------------------------------
# St_jacc_sphere — everything the shell RHS touches, staged on the device once.
#
# Untyped fields for the same reason as St_jacc2d: the concrete array type is
# whatever the active JACC backend hands back, and this is built before any of
# them is known. Read once per RHS call, never inside a kernel.
#---------------------------------------------------------------------------------
Base.@kwdef mutable struct St_jacc_sphere
    npoin::Int
    nelem::Int
    ngl::Int
    neq::Int

    # --- working precision (see jacc_working_float in rhs_jacc.jl) ----------
    TW
    mixed::Bool

    # --- mesh / metrics / basis, device-resident, read-only -----------------
    connijk
    coords
    Je
    dξdx
    dξdy
    dξdz
    dηdx
    dηdy
    dηdz
    dψ
    ω
    qe

    # --- point -> (element, i, j), CSR, for the race-free DSS ---------------
    p2e_ptr
    p2e_e
    p2e_i
    p2e_j

    # --- work ---------------------------------------------------------------
    q                    # npoin × (neq+1) device copy of the state
    q_stage              # host buffer in TW (mixed only)
    flux                 # nelem × ngl × ngl × 3neq   (equation-major)
    source               # nelem × ngl × ngl × neq
    gξ                   # nelem × ngl × ngl × neq    (viscous only)
    gη                   # nelem × ngl × ngl × neq    (viscous only)
    rhs_el               # nelem × ngl × ngl × neq
    RHS                  # npoin × neq
    RHS_stage            # host buffer in TW (mixed only)

    # --- run-time constants -------------------------------------------------
    μ                    # per-equation viscosity, device
    lvisc::Bool
    lpert::Bool
    PhysConst
end


#=================================================================================
                                  THE KERNELS
=================================================================================#

#--- node-local: the three flux components and the source ------------------------
function _jacc_sphere_volume!(iel, i, j, q, qe, coords, connijk,
                              flux, source, neq, PhysConst, lpert)

    @inbounds ip = connijk[iel, i, j]

    qip  = @view(q[ip, 1:neq])
    qeip = @view(qe[ip, :])

    # 3*neq values, equation-major: (F₁,G₁,H₁, F₂,G₂,H₂, …)
    @inbounds flux[iel, i, j, 1:3*neq] .= user_flux_gpu(qip, qeip, PhysConst, lpert)

    @inbounds source[iel, i, j, 1:neq] .= user_source_gpu(qip, qeip,
                                                          coords[1, ip], coords[2, ip],
                                                          coords[3, ip],
                                                          PhysConst, lpert)
    return nothing
end


#--- ν∇ₛq in the contravariant basis, with ωJ folded in --------------------------
#
# The serial _sphere_visc_el! applies ω[k]*ω[j]*Je[iel,k,j] to gξ inside the
# gather, and ω[i]*ω[k]*Je[iel,i,k] to gη. Both weights belong to the SOURCE node
# of their sum, so they can be folded in here and the gather below becomes a bare
# dψ-weighted sum — exactly the rearrangement the flat path makes.
#
function _jacc_sphere_visc_grad!(iel, i, j, q, connijk,
                                 dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                                 Je, dψ, ω, gξ, gη, μ, ngl, neq)

    T = eltype(gξ)

    @inbounds ωJac = ω[i]*ω[j]*Je[iel, i, j]

    @inbounds a1x = dξdx[iel, i, j]; a1y = dξdy[iel, i, j]; a1z = dξdz[iel, i, j]
    @inbounds a2x = dηdx[iel, i, j]; a2y = dηdy[iel, i, j]; a2z = dηdz[iel, i, j]

    @inbounds for ieq = 1:neq

        ν = μ[ieq]
        if ν <= zero(T)
            gξ[iel, i, j, ieq] = zero(T)
            gη[iel, i, j, ieq] = zero(T)
            continue
        end

        dqdξ = zero(T)
        dqdη = zero(T)
        for k = 1:ngl
            dqdξ += dψ[k, i]*q[connijk[iel, k, j], ieq]
            dqdη += dψ[k, j]*q[connijk[iel, i, k], ieq]
        end

        # ∇ₛq — tangential by construction, a¹ and a² being tangent vectors
        qx = dqdξ*a1x + dqdη*a2x
        qy = dqdξ*a1y + dqdη*a2y
        qz = dqdξ*a1z + dqdη*a2z

        gξ[iel, i, j, ieq] = ν*(qx*a1x + qy*a1y + qz*a1z)*ωJac
        gη[iel, i, j, ieq] = ν*(qx*a2x + qy*a2y + qz*a2z)*ωJac
    end
    return nothing
end


#--- the element residual at one node --------------------------------------------
#
#   rhs_el = -ωJ ( ∇ₛ·F - S )  -  ∫ ∇ₛψ·(ν∇ₛq)
#
# NOT divided by the mass matrix and NOT yet summed across the elements that share
# the node: _sphere_dss_scale! does both, on the host, after the cross-rank
# exchange — which is the order the serial path uses and the only order that is
# correct in parallel.
#
function _jacc_sphere_element_rhs!(iel, i, j, flux, source, gξ, gη,
                                   dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                                   Je, dψ, ω, rhs_el, ngl, neq, lvisc)

    T = eltype(rhs_el)

    @inbounds ωJac = ω[i]*ω[j]*Je[iel, i, j]

    @inbounds a1x = dξdx[iel, i, j]; a1y = dξdy[iel, i, j]; a1z = dξdz[iel, i, j]
    @inbounds a2x = dηdx[iel, i, j]; a2y = dηdy[iel, i, j]; a2z = dηdz[iel, i, j]

    @inbounds for ieq = 1:neq

        b = 3*(ieq - 1)          # equation-major flux base: F=b+1, G=b+2, H=b+3

        dFdξ = zero(T); dFdη = zero(T)
        dGdξ = zero(T); dGdη = zero(T)
        dHdξ = zero(T); dHdη = zero(T)

        for k = 1:ngl
            dψ_ki = dψ[k, i]
            dψ_kj = dψ[k, j]
            dFdξ += dψ_ki*flux[iel, k, j, b+1];  dFdη += dψ_kj*flux[iel, i, k, b+1]
            dGdξ += dψ_ki*flux[iel, k, j, b+2];  dGdη += dψ_kj*flux[iel, i, k, b+2]
            dHdξ += dψ_ki*flux[iel, k, j, b+3];  dHdη += dψ_kj*flux[iel, i, k, b+3]
        end

        dFdx = dFdξ*a1x + dFdη*a2x
        dGdy = dGdξ*a1y + dGdη*a2y
        dHdz = dHdξ*a1z + dHdη*a2z

        r = -ωJac*((dFdx + dGdy + dHdz) - source[iel, i, j, ieq])

        if lvisc
            t = zero(T)
            for k = 1:ngl
                t += dψ[i, k]*gξ[iel, k, j, ieq] + dψ[j, k]*gη[iel, i, k, ieq]
            end
            r -= t
        end

        rhs_el[iel, i, j, ieq] = r
    end
    return nothing
end


#--- direct stiffness summation, one thread per GLOBAL node ----------------------
#
# NO Minv here, unlike the flat path: on the shell the mass division happens in
# _sphere_dss_scale!, AFTER the cross-rank exchange. Applying it early would
# divide partial sums by the full mass on every shared node.
#
function _jacc_sphere_dss!(ip, rhs_el, RHS, p2e_ptr, p2e_e, p2e_i, p2e_j, neq)

    T = eltype(RHS)

    @inbounds lo = p2e_ptr[ip]
    @inbounds hi = p2e_ptr[ip+1] - 1

    @inbounds for ieq = 1:neq
        s = zero(T)
        for m = lo:hi
            s += rhs_el[p2e_e[m], p2e_i[m], p2e_j[m], ieq]
        end
        RHS[ip, ieq] = s
    end
    return nothing
end


#=================================================================================
                                    THE DRIVER
=================================================================================#

#---------------------------------------------------------------------------------
# sphere_rhs_jacc!(RHS, q, qe, mesh, metrics, sp, SVT)
#
# Drop-in for sphere_rhs!. `RHS` is the integrator's du (npoin × neqs+1 host
# matrix) and `q` is its state, so both stay on the host in Float64 and only the
# residual evaluation moves.
#---------------------------------------------------------------------------------
function sphere_rhs_jacc!(RHS, q, qe, mesh, metrics, sp, SVT)
    #
    # TYPED FUNCTION BARRIER, and it is not optional — the same lesson
    # _sphere_rhs_kernel! spells out a few hundred lines up in sphere_rhs.jl.
    #
    # sp.jacc is an untyped field (St_jacc_sphere is defined in this file, which
    # is included AFTER sphere_rhs.jl, so the struct cannot be named there). Read
    # without the assertion, `c` is ::Any, and then c.ngl and c.neq are ::Any too
    # — which makes `for k = 1:ngl` inside every kernel a dynamically-typed loop
    # and each iteration a dynamic dispatch. Measured on the Galewsky grid, 600
    # elements: 21.9 ms per RHS against 4.1 ms actually spent in the four kernels.
    # The assertion below is the whole difference.
    #
    return _sphere_rhs_jacc!(RHS, q, sp.jacc::St_jacc_sphere, mesh, metrics, sp)
end

function _sphere_rhs_jacc!(RHS, q, c::St_jacc_sphere, mesh, metrics, sp)

    npoin = c.npoin; nelem = c.nelem; ngl = c.ngl; neq = c.neq

    #--- state in: Float64 host → TW device ---------------------------------
    if c.mixed
        c.q_stage .= @view(q[:, 1:size(c.q_stage, 2)])
        copyto!(c.q, c.q_stage)
    else
        copyto!(c.q, @view(q[:, 1:size(c.q, 2)]))
    end

    JACC.parallel_for((nelem, ngl, ngl), _jacc_sphere_volume!,
                      c.q, c.qe, c.coords, c.connijk,
                      c.flux, c.source, neq, c.PhysConst, c.lpert)

    if c.lvisc
        JACC.parallel_for((nelem, ngl, ngl), _jacc_sphere_visc_grad!,
                          c.q, c.connijk,
                          c.dξdx, c.dξdy, c.dξdz, c.dηdx, c.dηdy, c.dηdz,
                          c.Je, c.dψ, c.ω, c.gξ, c.gη, c.μ, ngl, neq)
    end

    JACC.parallel_for((nelem, ngl, ngl), _jacc_sphere_element_rhs!,
                      c.flux, c.source, c.gξ, c.gη,
                      c.dξdx, c.dξdy, c.dξdz, c.dηdx, c.dηdy, c.dηdz,
                      c.Je, c.dψ, c.ω, c.rhs_el, ngl, neq, c.lvisc)

    JACC.parallel_for(npoin, _jacc_sphere_dss!,
                      c.rhs_el, c.RHS, c.p2e_ptr, c.p2e_e, c.p2e_i, c.p2e_j, neq)

    JACC.synchronize()

    #--- residual out: TW device → the integrator's Float64 du ---------------
    fill!(RHS, zero(eltype(RHS)))
    if c.mixed
        copyto!(c.RHS_stage, c.RHS)
        @view(RHS[:, 1:neq]) .= c.RHS_stage
    else
        copyto!(@view(RHS[:, 1:neq]), c.RHS)
    end

    # cross-rank DSS (a no-op on one rank) then M⁻¹ — the SERIAL function, so the
    # parallel semantics are identical by construction and not by argument.
    _sphere_dss_scale!(RHS, metrics.Minv, Int(mesh.npoin), Int(sp.neqs), sp.dss)

    return RHS
end


#---------------------------------------------------------------------------------
# build_jacc_cache_sphere(mesh, metrics, sp, inputs, SVT, T; rank)
#
# Stage the shell's mesh, metrics and reference state onto the device, once.
#---------------------------------------------------------------------------------
function build_jacc_cache_sphere(mesh, metrics, sp, qe, inputs, SVT, T; rank = 0)

    npoin = Int(mesh.npoin)
    nelem = Int(mesh.nelem)
    ngl   = Int(mesh.ngl)
    neq   = Int(sp.neqs)

    TW    = jacc_working_float(inputs; rank = rank)
    mixed = TW !== T

    get(inputs, :lfilter, false) == false || sp.lfilter == false ||
        @warn string("sphere_rhs_jacc.jl: :lfilter is on. The FILTER still runs on the ",
                     "host (sphere_filter!, once per step); only the RHS moves to JACC. ",
                     "That is correct but it costs a host round trip per step.")

    p2e_ptr, p2e_e, p2e_i, p2e_j = _jacc_build_p2e(mesh.connijk, npoin, nelem, ngl)

    vdim = sp.lvisc ? (nelem, ngl, ngl, neq) : (1, 1, 1, 1)

    return St_jacc_sphere(npoin = npoin, nelem = nelem, ngl = ngl, neq = neq,
                          TW = TW, mixed = mixed,

                          connijk = _jacc_stage(mesh.connijk, Int, (nelem, ngl, ngl)),
                          coords  = _jacc_stage(mesh.coords, TW, (NSD_MAX, npoin)),
                          Je      = _jacc_stage(metrics.Je,   TW, (nelem, ngl, ngl)),
                          dξdx    = _jacc_stage(metrics.dξdx, TW, (nelem, ngl, ngl)),
                          dξdy    = _jacc_stage(metrics.dξdy, TW, (nelem, ngl, ngl)),
                          dξdz    = _jacc_stage(metrics.dξdz, TW, (nelem, ngl, ngl)),
                          dηdx    = _jacc_stage(metrics.dηdx, TW, (nelem, ngl, ngl)),
                          dηdy    = _jacc_stage(metrics.dηdy, TW, (nelem, ngl, ngl)),
                          dηdz    = _jacc_stage(metrics.dηdz, TW, (nelem, ngl, ngl)),
                          dψ      = _jacc_stage(metrics.dψ,   TW, (ngl, ngl)),
                          ω       = _jacc_stage(metrics.ω,    TW, (ngl,)),
                          qe      = _jacc_stage(qe, TW, (npoin, size(qe, 2))),

                          p2e_ptr = _jacc_copy(p2e_ptr),
                          p2e_e   = _jacc_copy(p2e_e),
                          p2e_i   = _jacc_copy(p2e_i),
                          p2e_j   = _jacc_copy(p2e_j),

                          q         = _jacc_zeros(TW, npoin, neq+1),
                          q_stage   = mixed ? Base.zeros(TW, npoin, neq+1) : Base.zeros(TW, 0, 0),
                          flux      = _jacc_zeros(TW, nelem, ngl, ngl, 3*neq),
                          source    = _jacc_zeros(TW, nelem, ngl, ngl, neq),
                          gξ        = _jacc_zeros(TW, vdim...),
                          gη        = _jacc_zeros(TW, vdim...),
                          rhs_el    = _jacc_zeros(TW, nelem, ngl, ngl, neq),
                          RHS       = _jacc_zeros(TW, npoin, neq),
                          RHS_stage = mixed ? Base.zeros(TW, npoin, neq) : Base.zeros(TW, 0, 0),

                          μ         = _jacc_copy(TW.(sp.μ)),
                          lvisc     = sp.lvisc,
                          lpert     = SVT == PERT(),
                          PhysConst = TW === Float64 ? PHYS_CONST : PhysicalConst{TW}())
end


function jacc_sphere_banner(c::St_jacc_sphere; rank = 0)
    rank == 0 || return nothing
    nthr = Threads.nthreads()
    println(" # JACC shell RHS .................. backend = ", JACC.backend,
            JACC.backend == "threads" ?
                string("  (", nthr, " thread", nthr == 1 ? "" : "s", ")") : "",
            ",  arrays = ", nameof(typeof(c.RHS)))
    println(" #   npoin = ", c.npoin, ",  nelem = ", c.nelem, ",  ngl = ", c.ngl,
            ",  neqs = ", c.neq, ",  viscous = ", c.lvisc)
    println(" #   residual precision = ", c.TW,
            c.mixed ? "   (MIXED: state and MPI assembly stay in Float64)" : "")
    if JACC.backend == "threads" && nthr == 1
        @warn string("sphere_rhs_jacc.jl: :ljacc => true, but Julia has ONE thread, so the ",
                     "JACC path has no parallelism — and it is slower than the serial one, ",
                     "because JACC takes a different (unbatched) code path at ",
                     "nthreads() == 1. Start Julia with `-t <n>`.")
    end
    return nothing
end
