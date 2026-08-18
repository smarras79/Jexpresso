#---------------------------------------------------------------------------------
#
# rhs_jacc.jl — the 2D continuous-Galerkin SEM right-hand side, written ONCE
#               against JACC.jl and executed unchanged on
#
#                   threads (multi-core CPU) | CUDA | AMDGPU | oneAPI | Metal
#
# WHY A THIRD RHS AT ALL
# ----------------------
# The tree already carries two spellings of the same operator:
#
#   src/kernel/operators/rhs.jl      _build_rhs!            serial CPU, host arrays
#   src/kernel/operators/rhs_gpu.jl  _build_rhs_gpu_2D_v0!  KernelAbstractions
#
# The KA path works, but it is only reachable by setting `:backend => CUDABackend()`,
# which also flips TFloat/TInt to 32 bit globally and moves the WHOLE setup (mesh,
# metrics, bases, I/O) onto the device. That is a large, invasive switch, and every
# vendor needs its own KA backend package loaded before the deck is even read.
#
# JACC.jl inverts that: ONE kernel body, and the backend is a *build-time
# preference* (`JACC.set_backend("cuda")`, then restart Julia). The setup stays
# exactly where it is — on the host, in Float64 — and only the arrays the RHS
# actually reads are staged onto the device once, in `build_jacc_cache_2D`. On the
# `threads` backend `JACC.array` is the identity on a `Base.Array`, so the CPU run
# costs nothing but the staging copies and is a genuine multi-threaded CG solver.
#
# HOW IT IS TURNED ON
# -------------------
#   user_inputs.jl:  :ljacc => true          (with :backend left at its CPU default)
#   shell:           JULIA_NUM_THREADS=<n>   for the threads backend
#
# and, once per machine, for a GPU run:
#
#   julia> using JACC; JACC.set_backend("cuda")   # or "amdgpu" / "oneapi" / "metal"
#   (restart Julia — the backend is a Preferences setting, not a runtime switch)
#
# NO ATOMICS — AND THAT IS DELIBERATE
# -----------------------------------
# The KA kernels do the direct stiffness summation with `@atomic RHS[ip,ieq] -= …`
# from every (element, node) thread at once. That is correct but it is (a) not
# reproducible — the float sum order changes run to run, so two runs of the same
# deck do not agree bit for bit, and neither does a JACC run against the serial
# reference — and (b) it needs an atomic that every backend implements for the
# element type, which is exactly where portable code goes to die.
#
# Here the DSS is inverted instead: the element kernels write ONLY their own slot
# of `rhs_el[iel,i,j,ieq]` (no two threads share an output), and a second pass
# parallelised over GLOBAL NODES gathers the contributions through a precomputed
# point→(element,i,j) map:
#
#     RHS[ip,ieq] = Minv[ip] * Σ_{(iel,i,j) ∈ p2e[ip]} rhs_el[iel,i,j,ieq]
#
# Race-free by construction, deterministic to the last bit on every backend, and
# it costs one extra npoin×neqs pass plus an integer CSR map built once at setup.
# `_jacc_build_p2e` is that map.
#
# WHAT IS AND IS NOT COVERED
# --------------------------
#   yes   2D flat CG-SEM, `neqs` conservation laws, TOTAL() or PERT()
#   yes   the well-balanced volume term + user source
#   yes   constant artificial viscosity, :visc_model => AV()
#   yes   free-slip / Dirichlet boundaries through user_bc_dirichlet_gpu
#   yes   MPI (the cross-rank DSS is done on the host, exactly as the KA path does)
#   no    Laguerre semi-infinite elements, the CG filter, moisture, AMR
#   no    1D / 3D — see "Extending" at the bottom of this file
#   no    Neumann / surface-integral boundary conditions. The CPU path calls
#         apply_boundary_conditions_neumann! after the viscous term; there is no
#         counterpart here, so a case whose user_bc_neumann returns anything but
#         zero must NOT set :ljacc. This one cannot be detected from the deck —
#         the callback is only known by calling it — so it is a documented
#         restriction rather than a checked one.
#
# Anything in the "no" column raises at setup time rather than silently running a
# different equation.
#
# THE USER CALLBACKS ARE THE _gpu ONES
# ------------------------------------
# The kernels below call `user_flux_gpu`, `user_source_gpu`, `user_primitives_gpu`
# and `user_bc_dirichlet_gpu` from the case directory — the value-returning,
# allocation-free spellings the KA path already uses. A case that wants to run
# under JACC needs those four (three, if it is inviscid), not just the `!`-suffixed
# host versions.
#
#---------------------------------------------------------------------------------

using JACC


#---------------------------------------------------------------------------------
# Device-array staging.
#
# `JACC.array(x)` is `convert(JACC.array_type(), x)`. On the threads backend the
# array type IS `Base.Array`, so this returns `x` itself — the same object, not a
# copy. That aliasing is fine for the read-only inputs (mesh, metrics, basis) but
# NOT for the scratch buffers, which must be storage this file owns. `_jacc_copy`
# is therefore used for everything read-only and `_jacc_zeros` for scratch.
#---------------------------------------------------------------------------------
@inline _jacc_copy(a::AbstractArray) = JACC.array(collect(a))
@inline _jacc_zeros(T, dims...)      = JACC.zeros(T, dims...)


#---------------------------------------------------------------------------------
# St_jacc2d — everything the 2D RHS touches, staged on the device once.
#
# Deliberately UNTYPED fields: the concrete array type is whatever the active JACC
# backend hands back (Array / CuArray / ROCArray / …) and this struct is built
# before any of them is known. The fields are read exactly once per RHS call and
# passed straight into `JACC.parallel_for` as concrete arguments, so the dynamic
# lookup happens O(1) times per call and never inside a kernel.
#---------------------------------------------------------------------------------
Base.@kwdef mutable struct St_jacc2d
    # --- sizes -------------------------------------------------------------
    npoin::Int
    nelem::Int
    ngl::Int
    neq::Int
    nbnode::Int          # number of DISTINCT boundary nodes

    # --- mesh / metrics / basis (read-only on the device) -------------------
    connijk
    x
    y
    dξdx
    dξdy
    dηdx
    dηdy
    Je
    dψ
    ω
    Minv
    qe

    # --- point → (element, i, j) map, CSR, for the race-free DSS ------------
    p2e_ptr
    p2e_e
    p2e_i
    p2e_j

    # --- boundary node → (edge, local node) map, CSR, plus the normals ------
    bn_ip
    bn_ptr
    bn_edge
    bn_loc
    bdy_nx
    bdy_ny

    # --- work arrays --------------------------------------------------------
    u                    # npoin*neq, device copy of the ODE state
    uaux                 # npoin × (neq+1)
    flux                 # nelem × ngl × ngl × 2neq
    source               # nelem × ngl × ngl × neq
    uprim                # nelem × ngl × ngl × neq   (viscous only)
    gξ                   # nelem × ngl × ngl × neq   (viscous only)
    gη                   # nelem × ngl × ngl × neq   (viscous only)
    rhs_el               # nelem × ngl × ngl × neq
    qbdy                 # nbnode × neq
    RHS                  # npoin × neq
    RHS_host             # npoin × neq, host mirror handed to RHStoDU! / the MPI DSS

    # --- run-time constants -------------------------------------------------
    visc_coeff
    lvisc::Bool
    lsource::Bool
    lpert::Bool
    PhysConst
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
end


#---------------------------------------------------------------------------------
# _jacc_build_p2e(connijk, npoin, nelem, ngl)
#
# The inverse of `connijk`: for every global node ip, the list of (element, i, j)
# triples that map to it. Stored CSR-style so it can live on a GPU:
#
#     for m = ptr[ip] : ptr[ip+1]-1
#         iel = e[m];  i = ii[m];  j = jj[m]
#
# Built in linear time with a counting sort — the obvious `push!` into a
# Vector{Vector{Int}} is quadratic in memory traffic and allocates npoin vectors.
#---------------------------------------------------------------------------------
function _jacc_build_p2e(connijk, npoin::Int, nelem::Int, ngl::Int)

    ntot = nelem*ngl*ngl

    counts = zeros(Int, npoin)
    @inbounds for iel = 1:nelem, j = 1:ngl, i = 1:ngl
        counts[connijk[iel,i,j]] += 1
    end

    ptr = Vector{Int}(undef, npoin+1)
    ptr[1] = 1
    @inbounds for ip = 1:npoin
        ptr[ip+1] = ptr[ip] + counts[ip]
    end
    ptr[npoin+1] == ntot + 1 ||
        error(" # ERROR rhs_jacc.jl: connijk does not cover every node exactly once;",
              " got ", ptr[npoin+1]-1, " entries for ", ntot, " element nodes.")

    cursor = copy(ptr)
    pel = Vector{Int}(undef, ntot)     # not `pi`: that is Base.pi
    pix = Vector{Int}(undef, ntot)
    pjy = Vector{Int}(undef, ntot)
    @inbounds for iel = 1:nelem, j = 1:ngl, i = 1:ngl
        ip = connijk[iel,i,j]
        m  = cursor[ip]
        pel[m] = iel; pix[m] = i; pjy[m] = j
        cursor[ip] = m + 1
    end

    return ptr, pel, pix, pjy
end


#---------------------------------------------------------------------------------
# _jacc_build_bnodes(poin_in_bdy_edge, nedges_bdy, ngl)
#
# The same inversion for the boundary: the DISTINCT boundary nodes, each with the
# list of (edge, local node) slots it appears in.
#
# Why distinct nodes and not simply "one thread per (edge, node)", which is what
# apply_boundary_conditions_gpu! does: a corner node belongs to TWO edges, so the
# per-edge spelling has two threads writing the same u[ip] with two different
# normals — a race whose winner decides the corner's velocity. Here one thread
# owns the node and applies BOTH of its edges in a fixed order, which is
# deterministic and, at a free-slip corner, is also the physically right answer
# (both normal components get zeroed, not whichever one landed last).
#---------------------------------------------------------------------------------
function _jacc_build_bnodes(poin_in_bdy_edge, nedges_bdy::Int, ngl::Int)

    if nedges_bdy == 0
        return Int[], Int[1], Int[], Int[]
    end

    slots = Dict{Int, Vector{Tuple{Int,Int}}}()
    @inbounds for iedge = 1:nedges_bdy, k = 1:ngl
        ip = Int(poin_in_bdy_edge[iedge,k])
        ip > 0 || continue
        push!(get!(() -> Tuple{Int,Int}[], slots, ip), (iedge, k))
    end

    ips = sort!(collect(keys(slots)))          # sorted ⇒ the map is deterministic
    nb  = length(ips)

    ptr = Vector{Int}(undef, nb+1)
    ptr[1] = 1
    @inbounds for b = 1:nb
        ptr[b+1] = ptr[b] + length(slots[ips[b]])
    end

    edge = Vector{Int}(undef, ptr[nb+1]-1)
    loc  = Vector{Int}(undef, ptr[nb+1]-1)
    @inbounds for b = 1:nb
        m = ptr[b]
        for (iedge, k) in slots[ips[b]]
            edge[m] = iedge; loc[m] = k; m += 1
        end
    end

    return ips, ptr, edge, loc
end


#=================================================================================
                                  THE KERNELS

 Every one of them is a plain top-level function whose FIRST arguments are the
 loop indices JACC hands in, and whose remaining arguments are concrete arrays and
 isbits scalars. No closures, no captured globals, nothing that a GPU compiler
 cannot see through — that is the whole contract with `JACC.parallel_for`.
=================================================================================#

#--- u (flat, npoin*neq) → uaux (npoin × neq+1) -----------------------------------
function _jacc_u2uaux!(ip, u, uaux, npoin, neq)
    @inbounds for ieq = 1:neq
        uaux[ip, ieq] = u[(ieq-1)*npoin + ip]
    end
    return nothing
end


#--- Dirichlet / free-slip boundary values ----------------------------------------
#
# The sentinel protocol is the KA path's, kept verbatim so a case's
# user_bc_dirichlet_gpu does not need a second spelling: qbdy is pre-filled with
# 1234567 and the callback leaves that value in the slots it does not want to
# touch (SoliWaveIsland returns qbdy[1] for H — depth is free at a wall).
#
function _jacc_bc_2d!(ib, u, uaux, qe, x, y, t,
                      bn_ip, bn_ptr, bn_edge, bn_loc, bdy_nx, bdy_ny, qbdy,
                      neq, npoin, lpert)

    T = eltype(u)
    sentinel = T(1234567)

    @inbounds ip = bn_ip[ib]

    @inbounds for m = bn_ptr[ib]:(bn_ptr[ib+1]-1)

        iedge = bn_edge[m]
        k     = bn_loc[m]

        for ieq = 1:neq
            qbdy[ib, ieq] = sentinel
        end

        qbdy[ib, 1:neq] .= user_bc_dirichlet_gpu(@view(uaux[ip, :]), @view(qe[ip, :]),
                                                 x[ip], y[ip], t,
                                                 bdy_nx[iedge, k], bdy_ny[iedge, k],
                                                 @view(qbdy[ib, :]), lpert)

        for ieq = 1:neq
            qb = qbdy[ib, ieq]
            if qb != sentinel && qb != uaux[ip, ieq]
                uaux[ip, ieq]           = qb
                u[(ieq-1)*npoin + ip]   = qb
            end
        end
    end
    return nothing
end


#--- node-local quantities: F, G, S and (if viscous) the primitive variables -------
function _jacc_volume_2d!(iel, i, j, uaux, qe, x, y, connijk,
                          flux, source, uprim,
                          neq, lvisc, lsource, PhysConst, xmax, xmin, ymax, ymin, lpert)

    T = eltype(flux)

    @inbounds ip = connijk[iel, i, j]

    uip  = @view(uaux[ip, 1:neq])
    qeip = @view(qe[ip, 1:neq+1])

    @inbounds flux[iel, i, j, 1:2*neq] .= user_flux_gpu(uip, qeip, PhysConst, lpert)

    # `:lsource => false` means the deck has no source term at all — the CPU path
    # skips user_source! entirely in that case, so the JACC path must too, and not
    # merely add whatever a stale callback happens to return.
    if lsource
        @inbounds source[iel, i, j, 1:neq] .= user_source_gpu(uip, qeip, x[ip], y[ip],
                                                              PhysConst, xmax, xmin,
                                                              ymax, ymin, lpert)
    else
        @inbounds for ieq = 1:neq
            source[iel, i, j, ieq] = zero(T)
        end
    end

    if lvisc
        @inbounds uprim[iel, i, j, 1:neq] .= user_primitives_gpu(uip, @view(qe[ip, 1:neq]),
                                                                 lpert)
    end
    return nothing
end


#--- ν∇q in the contravariant basis, node by node ---------------------------------
#
# gξ = ωJ (a¹·ν∇q),  gη = ωJ (a²·ν∇q) with a¹ = (dξdx,dξdy), a² = (dηdx,dηdy).
# The ωJ factor is folded in HERE so the gather below is a bare dψ-weighted sum.
#
function _jacc_visc_grad_2d!(iel, i, j, uprim, dξdx, dξdy, dηdx, dηdy, Je, dψ, ω,
                             gξ, gη, visc_coeff, ngl, neq)

    T = eltype(uprim)

    @inbounds ωJac = ω[i]*ω[j]*Je[iel, i, j]

    @inbounds dξdx_ij = dξdx[iel, i, j]
    @inbounds dξdy_ij = dξdy[iel, i, j]
    @inbounds dηdx_ij = dηdx[iel, i, j]
    @inbounds dηdy_ij = dηdy[iel, i, j]

    @inbounds for ieq = 1:neq

        ν = visc_coeff[ieq]

        dqdξ = zero(T)
        dqdη = zero(T)
        for k = 1:ngl
            dqdξ += dψ[k, i]*uprim[iel, k, j, ieq]
            dqdη += dψ[k, j]*uprim[iel, i, k, ieq]
        end

        dqdx = ν*(dqdξ*dξdx_ij + dqdη*dηdx_ij)
        dqdy = ν*(dqdξ*dξdy_ij + dqdη*dηdy_ij)

        gξ[iel, i, j, ieq] = (dξdx_ij*dqdx + dξdy_ij*dqdy)*ωJac
        gη[iel, i, j, ieq] = (dηdx_ij*dqdx + dηdy_ij*dqdy)*ωJac
    end
    return nothing
end


#--- the element residual at one node ---------------------------------------------
#
# rhs_el = -ωJ ( ∇·F - S )  -  ∫∇ψ·(ν∇q)
#
# both terms in the form the gather pass expects: NOT yet divided by the mass
# matrix and NOT yet summed across the elements that share the node.
#
function _jacc_element_rhs_2d!(iel, i, j, flux, source, gξ, gη,
                               dξdx, dξdy, dηdx, dηdy, Je, dψ, ω,
                               rhs_el, ngl, neq, lvisc)

    T = eltype(rhs_el)

    @inbounds ωJac = ω[i]*ω[j]*Je[iel, i, j]

    @inbounds dξdx_ij = dξdx[iel, i, j]
    @inbounds dξdy_ij = dξdy[iel, i, j]
    @inbounds dηdx_ij = dηdx[iel, i, j]
    @inbounds dηdy_ij = dηdy[iel, i, j]

    @inbounds for ieq = 1:neq

        dFdξ = zero(T); dFdη = zero(T)
        dGdξ = zero(T); dGdη = zero(T)

        for k = 1:ngl
            dψ_ki = dψ[k, i]
            dψ_kj = dψ[k, j]
            dFdξ += dψ_ki*flux[iel, k, j, ieq]
            dFdη += dψ_kj*flux[iel, i, k, ieq]
            dGdξ += dψ_ki*flux[iel, k, j, neq + ieq]
            dGdη += dψ_kj*flux[iel, i, k, neq + ieq]
        end

        dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij
        dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij

        r = -ωJac*((dFdx + dGdy) - source[iel, i, j, ieq])

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


#--- direct stiffness summation + M⁻¹, one thread per GLOBAL node -----------------
function _jacc_dss_2d!(ip, rhs_el, RHS, Minv, p2e_ptr, p2e_e, p2e_i, p2e_j, neq)

    T = eltype(RHS)

    @inbounds mi = Minv[ip]
    @inbounds lo = p2e_ptr[ip]
    @inbounds hi = p2e_ptr[ip+1] - 1

    @inbounds for ieq = 1:neq
        s = zero(T)
        for m = lo:hi
            s += rhs_el[p2e_e[m], p2e_i[m], p2e_j[m], ieq]
        end
        RHS[ip, ieq] = s*mi
    end
    return nothing
end


#=================================================================================
                                    THE DRIVER
=================================================================================#

#---------------------------------------------------------------------------------
# jacc_rhs_2d!(c, u, t)
#
# Fills `c.RHS` (device, npoin × neq) with ∂q/∂t for the state `u`. `u` is the
# integrator's HOST vector; it is staged into `c.u` and — because the boundary
# kernel writes Dirichlet values back into it, exactly as the CPU and KA paths do
# — copied back out at the end.
#
# Everything after the staging is device-resident: six `parallel_for` launches and
# no host round trip.
#---------------------------------------------------------------------------------
function jacc_rhs_2d!(c::St_jacc2d, u, t)

    npoin = c.npoin; nelem = c.nelem; ngl = c.ngl; neq = c.neq
    T     = eltype(c.RHS)

    copyto!(c.u, u)

    JACC.parallel_for(npoin, _jacc_u2uaux!, c.u, c.uaux, npoin, neq)

    if c.nbnode > 0
        JACC.parallel_for(c.nbnode, _jacc_bc_2d!,
                          c.u, c.uaux, c.qe, c.x, c.y, T(t),
                          c.bn_ip, c.bn_ptr, c.bn_edge, c.bn_loc,
                          c.bdy_nx, c.bdy_ny, c.qbdy,
                          neq, npoin, c.lpert)
    end

    JACC.parallel_for((nelem, ngl, ngl), _jacc_volume_2d!,
                      c.uaux, c.qe, c.x, c.y, c.connijk,
                      c.flux, c.source, c.uprim,
                      neq, c.lvisc, c.lsource, c.PhysConst,
                      T(c.xmax), T(c.xmin), T(c.ymax), T(c.ymin), c.lpert)

    if c.lvisc
        JACC.parallel_for((nelem, ngl, ngl), _jacc_visc_grad_2d!,
                          c.uprim, c.dξdx, c.dξdy, c.dηdx, c.dηdy, c.Je, c.dψ, c.ω,
                          c.gξ, c.gη, c.visc_coeff, ngl, neq)
    end

    JACC.parallel_for((nelem, ngl, ngl), _jacc_element_rhs_2d!,
                      c.flux, c.source, c.gξ, c.gη,
                      c.dξdx, c.dξdy, c.dηdx, c.dηdy, c.Je, c.dψ, c.ω,
                      c.rhs_el, ngl, neq, c.lvisc)

    JACC.parallel_for(npoin, _jacc_dss_2d!,
                      c.rhs_el, c.RHS, c.Minv,
                      c.p2e_ptr, c.p2e_e, c.p2e_i, c.p2e_j, neq)

    JACC.synchronize()

    # The boundary kernel is the only writer of `u`; hand its Dirichlet values back
    # to the integrator's array, which is what both the CPU and the KA path do.
    copyto!(u, c.u)

    return c.RHS
end


#---------------------------------------------------------------------------------
# build_rhs_jacc_2D!(du, u, params, t)
#
# The `rhs!` entry point: run the device RHS, complete the cross-rank DSS on the
# host (a no-op on one rank), then unpack into `du` with the existing scalar-safe
# writer — `du` can be an OrdinaryDiffEq low-storage `ArrayFuse`, which supports
# only scalar setindex!, so it must not be written with a range or a broadcast.
#---------------------------------------------------------------------------------
function build_rhs_jacc_2D!(du, u, params, t)

    c = params.jacc

    jacc_rhs_2d!(c, u, t)

    copyto!(c.RHS_host, c.RHS)

    DSS_global_RHS!(@view(c.RHS_host[:,:]), params.g_dss_cache, params.neqs)

    RHStoDU!(du, @view(c.RHS_host[:,:]), params.neqs, params.mesh.npoin)

    return nothing
end


#=================================================================================
                                     SETUP
=================================================================================#

#---------------------------------------------------------------------------------
# jacc_check_inputs(inputs, sem, neqs)
#
# Fail at setup, loudly, on any configuration this file does not implement — a
# silently-ignored `:lfilter => true` is a different equation, not a slower run.
#---------------------------------------------------------------------------------
function jacc_check_inputs(inputs, mesh, VT)

    mesh.SD == NSD_2D() ||
        error(" # ERROR rhs_jacc.jl: :ljacc => true is implemented for 2D only.\n",
              " #   This case is ", mesh.SD, ". See \"Extending\" at the bottom of",
              " src/kernel/operators/rhs_jacc.jl.")

    inputs[:backend] == CPU() ||
        error(" # ERROR rhs_jacc.jl: :ljacc => true wants :backend => CPU() (the default).\n",
              " #   JACC owns the device here; :backend selects the KernelAbstractions\n",
              " #   path instead, and the two must not both be on. Pick the JACC device\n",
              " #   with `using JACC; JACC.set_backend(\"cuda\")` and restart Julia.")

    mesh.lLaguerre == false &&
        get(inputs, :llaguerre_1d_right, false) == false &&
        get(inputs, :llaguerre_1d_left,  false) == false ||
        error(" # ERROR rhs_jacc.jl: the JACC path has no Laguerre semi-infinite elements.")

    get(inputs, :lspherical_shell, false) == false ||
        error(" # ERROR rhs_jacc.jl: :lspherical_shell => true does not go through rhs! at all —\n",
              " #   the shell has its own RHS (sphere_rhs!) and its own time loop. Porting it\n",
              " #   is the next step; see \"Extending\" at the bottom of this file.")

    get(inputs, :ladapt, false) == false && get(inputs, :lamr, false) == false ||
        error(" # ERROR rhs_jacc.jl: the JACC cache is staged onto the device ONCE, at setup.\n",
              " #   An adaptive run (:ladapt / :lamr) remeshes underneath it, so the cache\n",
              " #   would go stale silently. Not supported.")

    get(inputs, :lfilter, false) == false ||
        error(" # ERROR rhs_jacc.jl: the JACC path has no CG filter yet (:lfilter => true).")

    get(inputs, :lmoist, false) == false ||
        error(" # ERROR rhs_jacc.jl: the JACC path has no microphysics yet (:lmoist => true).")

    if get(inputs, :lvisc, false) == true
        VT == AV() ||
            error(" # ERROR rhs_jacc.jl: the JACC path implements :visc_model => AV() only;\n",
                  " #   this case asks for ", VT, ".")
    end

    return nothing
end


#---------------------------------------------------------------------------------
# build_jacc_cache_2D(sem, qp, inputs, T)
#
# Stage the mesh, the metrics, the basis and the reference state onto the device
# and allocate the work arrays. Called ONCE, from params_setup.
#---------------------------------------------------------------------------------
function build_jacc_cache_2D(sem, qp, inputs, T)

    mesh    = sem.mesh
    metrics = sem.metrics
    basis   = sem.basis
    ω       = sem.ω

    npoin = Int(mesh.npoin)
    nelem = Int(mesh.nelem)
    ngl   = Int(mesh.ngl)
    neq   = Int(qp.neqs)

    jacc_check_inputs(inputs, mesh, get(inputs, :visc_model, AV()))

    lvisc = get(inputs, :lvisc, false) == true
    lpert = sem.SOL_VARS_TYPE == PERT()

    #--- the two CSR maps that replace the atomics ---------------------------
    p2e_ptr, p2e_e, p2e_i, p2e_j = _jacc_build_p2e(mesh.connijk, npoin, nelem, ngl)

    nedges_bdy = Int(mesh.nedges_bdy)
    bn_ip, bn_ptr, bn_edge, bn_loc = _jacc_build_bnodes(mesh.poin_in_bdy_edge,
                                                        nedges_bdy, ngl)
    nbnode = length(bn_ip)

    #--- per-equation viscosity ---------------------------------------------
    coeffs = zeros(T, neq)
    if lvisc
        μ = inputs[:μ]
        for ieq = 1:neq
            coeffs[ieq] = T(μ isa Number ? μ : μ[min(ieq, length(μ))])
        end
    end

    #--- viscous scratch is 1×1×1×1 when the run is inviscid ----------------
    vdim = lvisc ? (nelem, ngl, ngl, neq) : (1, 1, 1, 1)

    return St_jacc2d(npoin = npoin, nelem = nelem, ngl = ngl, neq = neq,
                     nbnode = nbnode,

                     connijk = _jacc_copy(Int.(mesh.connijk)),
                     x       = _jacc_copy(T.(mesh.x)),
                     y       = _jacc_copy(T.(mesh.y)),
                     dξdx    = _jacc_copy(T.(metrics.dξdx)),
                     dξdy    = _jacc_copy(T.(metrics.dξdy)),
                     dηdx    = _jacc_copy(T.(metrics.dηdx)),
                     dηdy    = _jacc_copy(T.(metrics.dηdy)),
                     Je      = _jacc_copy(T.(metrics.Je)),
                     dψ      = _jacc_copy(T.(basis.dψ)),
                     ω       = _jacc_copy(T.(ω)),
                     Minv    = _jacc_copy(T.(sem.matrix.Minv)),
                     qe      = _jacc_copy(T.(qp.qe)),

                     p2e_ptr = _jacc_copy(p2e_ptr),
                     p2e_e   = _jacc_copy(p2e_e),
                     p2e_i   = _jacc_copy(p2e_i),
                     p2e_j   = _jacc_copy(p2e_j),

                     bn_ip   = _jacc_copy(bn_ip),
                     bn_ptr  = _jacc_copy(bn_ptr),
                     bn_edge = _jacc_copy(bn_edge),
                     bn_loc  = _jacc_copy(bn_loc),
                     bdy_nx  = _jacc_copy(nedges_bdy > 0 ? T.(metrics.nx) : zeros(T,1,1)),
                     bdy_ny  = _jacc_copy(nedges_bdy > 0 ? T.(metrics.ny) : zeros(T,1,1)),

                     u        = _jacc_zeros(T, npoin*neq),
                     uaux     = _jacc_zeros(T, npoin, neq+1),
                     flux     = _jacc_zeros(T, nelem, ngl, ngl, 2*neq),
                     source   = _jacc_zeros(T, nelem, ngl, ngl, neq),
                     uprim    = _jacc_zeros(T, vdim...),
                     gξ       = _jacc_zeros(T, vdim...),
                     gη       = _jacc_zeros(T, vdim...),
                     rhs_el   = _jacc_zeros(T, nelem, ngl, ngl, neq),
                     qbdy     = _jacc_zeros(T, max(nbnode,1), neq),
                     RHS      = _jacc_zeros(T, npoin, neq),
                     RHS_host = zeros(T, npoin, neq),

                     visc_coeff = _jacc_copy(coeffs),
                     lvisc = lvisc,
                     lsource = get(inputs, :lsource, false) == true,
                     lpert = lpert,
                     PhysConst = PHYS_CONST,
                     xmin = Float64(mesh.xmin), xmax = Float64(mesh.xmax),
                     ymin = Float64(mesh.ymin), ymax = Float64(mesh.ymax))
end


#---------------------------------------------------------------------------------
# jacc_banner(c) — one line at startup saying which device is actually in use.
#
# JACC's backend is a Preferences value baked in at load time, so this is the only
# place a user finds out whether the `JACC.set_backend("cuda")` they ran last week
# in a different project directory is the one this run picked up.
#---------------------------------------------------------------------------------
function jacc_banner(c::St_jacc2d; rank = 0)
    rank == 0 || return nothing
    nthr = Threads.nthreads()
    println(" # JACC RHS ........................ backend = ", JACC.backend,
            JACC.backend == "threads" ? string("  (", nthr, " thread", nthr == 1 ? "" : "s", ")") : "",
            ",  arrays = ", nameof(typeof(c.RHS)))
    println(" #   npoin = ", c.npoin, ",  nelem = ", c.nelem, ",  ngl = ", c.ngl,
            ",  neqs = ", c.neq, ",  viscous = ", c.lvisc)
    return nothing
end


#=================================================================================
 EXTENDING THIS FILE

 1D / 3D
   The only shape-dependent pieces are `_jacc_build_p2e` (loop nest over the
   element's local indices), `_jacc_volume_2d!` (which coordinates go to the user
   callbacks), `_jacc_element_rhs_2d!` (how many contravariant directions) and the
   `parallel_for` launch dims. Everything else — the CSR gather, the boundary node
   map, the cache, the driver — is dimension-agnostic as written.

 The spherical shell (problems/ShallowWater/SWsphere)
   The shell does NOT go through `rhs!`; it has its own RHS
   (`sphere_rhs!`, src/kernel/operators/sphere_rhs.jl) and its own time loop
   (src/kernel/solvers/sphere_time_loop.jl). Porting it is the same three kernels
   as here with one change of substance: the manifold metric is 3×2, so a¹ and a²
   each carry a z component and the divergence is dFdx + dGdy + dHdz with a THIRD
   flux component H. `_jacc_element_rhs_2d!` grows one term; the gather and the
   boundary map are unchanged (a closed shell has no boundary at all, so nbnode is
   zero and the BC kernel never launches).

 The CG filter, moisture, Laguerre
   Not started. `jacc_check_inputs` refuses them by name rather than running a
   different equation quietly.
=================================================================================#
