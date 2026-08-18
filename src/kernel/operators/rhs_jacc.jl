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
# Race-free by construction, deterministic to the last bit on every backend and at
# every thread count, and it costs one extra npoin×neqs pass plus an integer CSR
# map built once at setup. `_jacc_build_p2e` is that map.
#
# To be exact about what that determinism buys: two runs of THIS path agree
# bit for bit. Against the CPU path the agreement is to ROUND-OFF, not to the bit,
# precisely because the summation order differs — a gather per node adds the same
# contributions in a different sequence than a scatter per element, and floating
# point addition is not associative. Measured on SoliWaveIsland, five steps of
# SSPRK54: 4.4e-16 relative on H, 7.7e-16 on Hu.
#
# WORKING PRECISION, AND WHY IT IS A SEPARATE KNOB
# ------------------------------------------------
# Jexpresso carries its state in Float64. Some devices cannot: Apple GPUs have no
# double precision at all — Metal.jl's array constructor refuses Float64 outright
# ("Metal does not support Float64 values"), and JACC agrees with it, declaring
# `default_float(::MetalBackend) = Float32`.
#
# So the device arrays here carry their own float type, `TW`, chosen by
# `:jacc_float` (default: whatever the active backend can do). When TW is Float32
# the scheme is MIXED precision, not single precision:
#
#     state u, du, the MPI assembly, the time integrator   Float64   (untouched)
#     the residual evaluation — fluxes, source, viscosity  Float32   (on device)
#
# That is the useful arrangement, and it is not the same as running the solver in
# single precision. The state never loses digits to repeated rounding: it is read
# down to Float32 once per call and the Float32 residual is added back into a
# Float64 state, so the error is O(eps32) per step in the TENDENCY rather than
# O(eps32) accumulating in the solution. What it costs is the accuracy of the
# residual itself: on this deck the solitary wave is η ≈ 6e-4 riding on
# h0 = 0.32, a relative perturbation of 2e-3, so with eps32 ≈ 1.2e-7 about four
# significant digits of the signal survive. Fine for a demo or a throughput
# measurement, and NOT something to trust for a vorticity-sensitive run.
#
# The lake-at-rest equilibrium does survive the conversion exactly: H and He are
# equal before the cast, so they are equal after it, and the well-balanced
# flux/source split still cancels term by term. That is checked, in Float32, by
# test/test_jacc_rhs.jl.
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

#
# Bring the vendor package for the ACTIVE backend into scope.
#
# JACC.set_backend("cuda") adds CUDA.jl to the project's [deps] and records the
# choice in LocalPreferences.toml, but it does not import it — and the CUDA
# array type, the CUDA parallel_for and everything else backend-specific live in
# a JACC package EXTENSION that only loads once CUDA itself is imported. Without
# this line a GPU run gets a MethodError on `JACC.array_type()` at cache-build
# time instead of a GPU.
#
# On the default "threads" backend the macro expands to a single @info line and
# imports nothing. It has to sit at top level (a macro expands where it is
# written), and it is recompiled automatically when the preference changes:
# JACC reads the backend with @load_preference, so Preferences.jl invalidates
# JACC's cache — and transitively Jexpresso's — the moment set_backend runs.
#
JACC.@init_backend


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
# _jacc_stage(a, T, dims) — stage `a` on the device as element type T and EXACTLY
# rank `length(dims)`.
#
# WHY THE RESHAPE IS NOT COSMETIC. Jexpresso's 2D metric arrays carry trailing
# singleton dimensions: allocate_metrics (metric_terms.jl) gives dξdx, dηdy, Je …
# the shape (nelem, ngl, ngl, 1) and the boundary normals nx, ny the shape
# (nedges_bdy, ngl, 1). The rest of the tree indexes them with one index fewer —
# `Je[iel,i,j]`, `nx[iedge,k]` — and gets away with it because Julia lets a
# Base.Array drop TRAILING SINGLETON dimensions.
#
# That rule does not survive the trip into a kernel. JACC's threads backend runs
# the body under Polyester.@batch, which hands the arguments in as
# StrideArraysCore.PtrArray, and a PtrArray demands exactly as many indices as it
# has dimensions: `Je[iel,i,j]` on a 4-D array is a BoundsError, not a read. On a
# GPU backend the device array types are just as strict.
#
# This is not hypothetical. It is what the first full-solver run of this file did,
# and only with more than one thread — at nthreads() == 1 JACC takes a plain loop,
# the arrays stay Base.Arrays, and the same code is fine. A bug that appears only
# under threading is worth removing at the source, so the cache stages every array
# at the rank its kernel actually indexes. The copies are ours, so it is free.
#---------------------------------------------------------------------------------
function _jacc_stage(a, ::Type{T}, dims::NTuple{N,Int}) where {T,N}
    length(a) == prod(dims) ||
        error(" # ERROR rhs_jacc.jl: cannot stage an array of size ", size(a),
              " as ", dims, " — element counts differ (", length(a),
              " vs ", prod(dims), ").")
    return _jacc_copy(T.(reshape(collect(a), dims)))
end


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
    coords               # nsd × npoin — the canonical coordinate storage
    x                    # kept: user_source_gpu still takes x, y as scalars
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

    # --- working precision --------------------------------------------------
    TW                   # element type of every device float array
    mixed::Bool          # TW !== Float64, so the host↔device hops must convert
    u_stage              # npoin*neq host buffer in TW   (mixed only)
    bout                 # nbnode × neq device: the corrected boundary values
    bout_host            # nbnode × neq host mirror of the above
    bn_ip_host           # nbnode host copy of bn_ip, to scatter bout into u
    RHS_stage            # npoin × neq host buffer in TW (mixed only)

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



#---------------------------------------------------------------------------------
# _jacc_almost_equal(a, b) — Kopriva's algorithm 139, the SAME test
# build_custom_bcs_dirichlet! uses to decide whether a Dirichlet value is worth
# writing (AlmostEqual, src/kernel/infrastructure/Kopriva_functions.jl).
#
# Note the tolerance: ε = 1e-6, ABSOLUTE for small values. It is enormous next to
# a physical momentum, and it is not an accident of this port — it is what the CPU
# path does, so any boundary correction smaller than ~1e-6 is silently declined
# there. Reproducing it is the difference between a JACC run that tracks the
# reference and one that quietly enforces a slightly different wall.
#
# This was measured, not guessed: with an exact `!=` here, a corner node of the
# SoliWaveIsland basin carried Hu = 8.95e-7 on the CPU path and exactly 0 on this
# one, and the two runs parted company at the fourth time step.
#---------------------------------------------------------------------------------
@inline function _jacc_almost_equal(a::T, b::T) where {T}
    ε = T(1.0e-6)
    if a == zero(T) || b == zero(T) || a <= ε || b <= ε
        return abs(a - b) <= 2ε
    else
        return abs(a - b) <= ε*abs(a) && abs(a - b) <= ε*abs(b)
    end
end


#--- Dirichlet / free-slip boundary values ----------------------------------------
#
# The sentinel protocol is the KA path's, so a case's user_bc_dirichlet_gpu needs
# no second spelling: qbdy is pre-filled with a sentinel and the callback leaves it
# in the slots it does not want to touch (SoliWaveIsland returns qbdy[1] for H —
# the depth is free at a wall).
#
# THE BOUNDARY VALUES DO REACH `u`, and that is not optional.
#
# build_custom_bcs_dirichlet! (BCs.jl) opens with
#
#   # WARNING: Notice that the b.c. are applied to uaux[:,:] and NOT u[:]!
#
# and then closes, as its very last statement, with `uaux2u!(u, uaux, neqs, npoin)`.
# The comment is stale: the corrected field IS pushed back into the integrator's
# state on every RHS call. Measured on this deck: the CPU path changes `u` by
# 1.06e-3 during the first RK stage, all of it at wall nodes.
#
# So the kernel corrects `uaux` and the driver mirrors it into `u` afterwards
# (_jacc_uaux2u!), exactly as the CPU does. Skipping it is not a rounding-level
# difference — with the write left out the two paths part company by 1e-3 inside
# the first time step, because the wall is then enforced in one trajectory and not
# in the other.
#
function _jacc_bc_2d!(ib, uaux, qe, coords, t,
                      bn_ip, bn_ptr, bn_edge, bn_loc, bdy_nx, bdy_ny, qbdy, bout,
                      neq, lpert)

    T = eltype(uaux)
    sentinel = T(1234567)

    @inbounds ip = bn_ip[ib]

    @inbounds for m = bn_ptr[ib]:(bn_ptr[ib+1]-1)

        iedge = bn_edge[m]
        k     = bn_loc[m]

        for ieq = 1:neq
            qbdy[ib, ieq] = sentinel
        end

        qbdy[ib, 1:neq] .= user_bc_dirichlet_gpu(@view(uaux[ip, :]), @view(qe[ip, :]),
                                                 @view(coords[:, ip]), t,
                                                 bdy_nx[iedge, k], bdy_ny[iedge, k],
                                                 @view(qbdy[ib, :]), lpert)

        for ieq = 1:neq
            qb = qbdy[ib, ieq]
            if !_jacc_almost_equal(qb, sentinel) && !_jacc_almost_equal(qb, uaux[ip, ieq])
                uaux[ip, ieq] = qb
            end
        end
    end

    # The corrected values, packed nbnode-wide so the driver can hand them back
    # to the integrator with a transfer the size of the BOUNDARY rather than the
    # size of the field. Away from the boundary uaux2u! is the identity, so this
    # carries exactly the information the CPU's whole-field write carries.
    @inbounds for ieq = 1:neq
        bout[ib, ieq] = uaux[ip, ieq]
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
# integrator's HOST vector; it is staged into `c.u` on the way in, and the
# boundary-corrected field is handed back on the way out — the same contract the
# CPU path has (see _jacc_bc_2d!).
#
# Everything in between is device-resident: five to seven `parallel_for` launches
# and no host round trip.
#---------------------------------------------------------------------------------
function jacc_rhs_2d!(c::St_jacc2d, u, t)

    npoin = c.npoin; nelem = c.nelem; ngl = c.ngl; neq = c.neq
    T     = c.TW

    #--- state in: Float64 host → TW device ---------------------------------
    # The `mixed` branch is the whole cost of reduced precision on the way in:
    # one host-side cast into a preallocated buffer, then the usual transfer. In
    # the Float64 case there is no cast and no extra buffer, so the common path
    # is exactly what it was.
    if c.mixed
        c.u_stage .= u
        copyto!(c.u, c.u_stage)
    else
        copyto!(c.u, u)
    end

    JACC.parallel_for(npoin, _jacc_u2uaux!, c.u, c.uaux, npoin, neq)

    if c.nbnode > 0
        JACC.parallel_for(c.nbnode, _jacc_bc_2d!,
                          c.uaux, c.qe, c.coords, T(t),
                          c.bn_ip, c.bn_ptr, c.bn_edge, c.bn_loc,
                          c.bdy_nx, c.bdy_ny, c.qbdy, c.bout,
                          neq, c.lpert)
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

    #--- boundary values back into the integrator's state -------------------
    #
    # uaux2u! at the tail of the CPU's Dirichlet builder writes the WHOLE field
    # back. Away from the boundary that writes back exactly what it read, so
    # scattering the nbnode corrected rows is the same operation — and it moves
    # nbnode×neq numbers instead of npoin×neq (≈440 vs ≈37000 on this deck).
    #
    # It also keeps mixed precision honest. A whole-field write-back would round
    # every node of a Float64 state through Float32 on every RHS call, which is
    # the one thing the mixed scheme is supposed to avoid; only the boundary
    # nodes, whose values the device actually computed, are demoted.
    #
    if c.nbnode > 0
        copyto!(c.bout_host, c.bout)
        @inbounds for ieq = 1:neq
            off = (ieq-1)*npoin
            for ib = 1:c.nbnode
                u[off + c.bn_ip_host[ib]] = c.bout_host[ib, ieq]
            end
        end
    end

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

    #--- RHS out: TW device → Float64 host -----------------------------------
    if c.mixed
        copyto!(c.RHS_stage, c.RHS)
        c.RHS_host .= c.RHS_stage
    else
        copyto!(c.RHS_host, c.RHS)
    end

    #
    # params.uaux is the CPU path's scratch, and things downstream of rhs! read it
    # as if rhs! had just refreshed it — _build_rhs! fills it from u on every call,
    # so diagnostics, the LES statistics and some of the writers assume it holds
    # the current state plus the boundary values. Without this they silently see
    # whatever the last CPU-path call left behind (the initial condition, on a run
    # that never took the CPU path at all).
    #
    # Rebuilt from `u` rather than copied off the device: jacc_rhs_2d! has just
    # scattered the boundary values into `u`, so u2uaux! reproduces exactly the
    # field the device kernels worked on — for the price of one host pass and no
    # transfer at all.
    #
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)

    DSS_global_RHS!(@view(c.RHS_host[:,:]), params.g_dss_cache, params.neqs)

    RHStoDU!(du, @view(c.RHS_host[:,:]), params.neqs, params.mesh.npoin)

    return nothing
end


#=================================================================================
                                     SETUP
=================================================================================#

#---------------------------------------------------------------------------------
# jacc_working_float(inputs; rank) — which float type the DEVICE arrays carry.
#
#   :jacc_float => :auto      (default) Float64 where the backend can do it, and
#                             the backend's own default where it cannot.
#   :jacc_float => Float32    ask for it explicitly, no warning.
#   :jacc_float => Float64    demand it; an error if the backend has no Float64.
#
# `JACC.default_float()` is the backend's own declaration — Float64 for threads,
# CUDA, AMDGPU and oneAPI, Float32 for Metal, whose hardware has no double
# precision at all.
#
# The :auto fallback WARNS rather than falling back silently. Quietly halving the
# precision of a geophysical solve because of the machine it happened to be
# launched on is not a convenience.
#---------------------------------------------------------------------------------
function jacc_working_float(inputs; rank = 0)

    backend_float = JACC.default_float()
    want = get(inputs, :jacc_float, :auto)

    if want === :auto || want == "auto"
        if backend_float !== Float64 && rank == 0
            @warn string("rhs_jacc.jl: the active JACC backend (\"", JACC.backend,
                         "\") has no Float64, so the RHS will be evaluated in ",
                         backend_float, ".\n",
                         "The state, the time integrator and the MPI assembly stay in ",
                         "Float64 — this is a MIXED-precision residual, not a ",
                         "single-precision solver — but the tendency now carries only ",
                         "~", round(Int, -log10(eps(backend_float))), " significant digits.\n",
                         "Set :jacc_float => Float32 in the deck to accept this ",
                         "silently, or run on a backend with double precision.")
        end
        return backend_float
    end

    TW = want isa Type ? want :
         (want in (:f32, "f32", :Float32, "Float32") ? Float32 :
          want in (:f64, "f64", :Float64, "Float64") ? Float64 : nothing)

    TW in (Float32, Float64) ||
        error(" # ERROR rhs_jacc.jl: :jacc_float must be :auto, Float32 or Float64;\n",
              " #   got ", repr(want), ".")

    TW === Float64 && backend_float !== Float64 &&
        error(" # ERROR rhs_jacc.jl: the deck asks for :jacc_float => Float64, but the\n",
              " #   active JACC backend (\"", JACC.backend, "\") has no double precision",
              " — JACC\n #   itself reports default_float() = ", backend_float,
              ", and Metal.jl refuses a Float64\n",
              " #   array outright. Use :jacc_float => Float32 (mixed precision: the\n",
              " #   state stays Float64, only the residual is evaluated in single), or\n",
              " #   pick a backend that has Float64.")

    return TW
end


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
function build_jacc_cache_2D(sem, qp, inputs, T; rank = 0)

    #
    # T is the HOST float type (TFloat, i.e. Float64). TW is what the device
    # arrays carry — the same thing on every backend that has double precision,
    # Float32 on Metal. See jacc_working_float.
    #
    TW    = jacc_working_float(inputs; rank = rank)
    mixed = TW !== T

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
    #
    # Built from host copies reshaped to the rank the loops below index, for the
    # same reason _jacc_stage exists: mesh.connijk and mesh.poin_in_bdy_edge may
    # arrive with a trailing singleton dimension.
    #
    conn_h = reshape(collect(mesh.connijk), nelem, ngl, ngl)
    p2e_ptr, p2e_e, p2e_i, p2e_j = _jacc_build_p2e(conn_h, npoin, nelem, ngl)

    nedges_bdy = Int(mesh.nedges_bdy)
    bdy_h  = reshape(collect(mesh.poin_in_bdy_edge), nedges_bdy, ngl)
    bn_ip, bn_ptr, bn_edge, bn_loc = _jacc_build_bnodes(bdy_h, nedges_bdy, ngl)
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

                     TW      = TW,
                     mixed   = mixed,
                     u_stage    = mixed ? Base.zeros(TW, npoin*neq) : Base.zeros(TW, 0),
                     bout       = _jacc_zeros(TW, max(nbnode,1), neq),
                     bout_host  = Base.zeros(TW, max(nbnode,1), neq),
                     bn_ip_host = bn_ip,
                     RHS_stage  = mixed ? Base.zeros(TW, npoin, neq) : Base.zeros(TW, 0, 0),

                     connijk = _jacc_stage(mesh.connijk, Int, (nelem, ngl, ngl)),
                     coords  = _jacc_stage(mesh.coords, TW, (Int(size(mesh.coords,1)), npoin)),
                     x       = _jacc_stage(mesh.x,    TW, (npoin,)),
                     y       = _jacc_stage(mesh.y,    TW, (npoin,)),
                     dξdx    = _jacc_stage(metrics.dξdx, TW, (nelem, ngl, ngl)),
                     dξdy    = _jacc_stage(metrics.dξdy, TW, (nelem, ngl, ngl)),
                     dηdx    = _jacc_stage(metrics.dηdx, TW, (nelem, ngl, ngl)),
                     dηdy    = _jacc_stage(metrics.dηdy, TW, (nelem, ngl, ngl)),
                     Je      = _jacc_stage(metrics.Je,   TW, (nelem, ngl, ngl)),
                     dψ      = _jacc_stage(basis.dψ,     TW, (ngl, ngl)),
                     ω       = _jacc_stage(ω,            TW, (ngl,)),
                     Minv    = _jacc_stage(sem.matrix.Minv, TW, (npoin,)),
                     qe      = _jacc_stage(qp.qe, TW, (npoin, neq+1)),

                     p2e_ptr = _jacc_copy(p2e_ptr),
                     p2e_e   = _jacc_copy(p2e_e),
                     p2e_i   = _jacc_copy(p2e_i),
                     p2e_j   = _jacc_copy(p2e_j),

                     bn_ip   = _jacc_copy(bn_ip),
                     bn_ptr  = _jacc_copy(bn_ptr),
                     bn_edge = _jacc_copy(bn_edge),
                     bn_loc  = _jacc_copy(bn_loc),
                     bdy_nx  = nedges_bdy > 0 ?
                               _jacc_stage(metrics.nx, TW, (nedges_bdy, ngl)) :
                               _jacc_zeros(TW, 1, 1),
                     bdy_ny  = nedges_bdy > 0 ?
                               _jacc_stage(metrics.ny, TW, (nedges_bdy, ngl)) :
                               _jacc_zeros(TW, 1, 1),

                     u        = _jacc_zeros(TW, npoin*neq),
                     uaux     = _jacc_zeros(TW, npoin, neq+1),
                     flux     = _jacc_zeros(TW, nelem, ngl, ngl, 2*neq),
                     source   = _jacc_zeros(TW, nelem, ngl, ngl, neq),
                     uprim    = _jacc_zeros(TW, vdim...),
                     gξ       = _jacc_zeros(TW, vdim...),
                     gη       = _jacc_zeros(TW, vdim...),
                     rhs_el   = _jacc_zeros(TW, nelem, ngl, ngl, neq),
                     qbdy     = _jacc_zeros(TW, max(nbnode,1), neq),
                     RHS      = _jacc_zeros(TW, npoin, neq),
                     RHS_host = Base.zeros(T, npoin, neq),

                     visc_coeff = _jacc_copy(TW.(coeffs)),
                     lvisc = lvisc,
                     lsource = get(inputs, :lsource, false) == true,
                     lpert = lpert,
                     PhysConst = TW === Float64 ? PHYS_CONST : PhysicalConst{TW}(),
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
    println(" #   residual precision = ", c.TW,
            c.mixed ? "   (MIXED: state, integrator and MPI assembly stay in Float64)" :
                      "")

    #
    # The single-thread trap. `:ljacc => true` says WHICH kernels run; it does not
    # give Julia any threads — Julia starts with one unless -t / JULIA_NUM_THREADS
    # says otherwise, and nothing in Jexpresso changes that.
    #
    # It is worth a warning and not just a doc line because the cost is not the
    # obvious one. JACC's threads backend runs a plain loop at nthreads() == 1 and
    # Polyester.@batch above it, and the two are not the same code. Measured on
    # this RHS, 32x16 elements at nop=5: 51.5 RHS/s on one thread against 397.6 on
    # two — 7.7x, which is not something two threads can do. So a one-thread JACC
    # run is slower than leaving :ljacc off altogether.
    #
    if JACC.backend == "threads" && nthr == 1
        @warn string("rhs_jacc.jl: :ljacc => true, but Julia is running with ONE thread, ",
                     "so the JACC path has no parallelism at all — and it is slower than ",
                     "the serial CPU path, because JACC takes a different (unbatched) ",
                     "code path at nthreads() == 1.\n",
                     "Start Julia with `-t <n>` or JULIA_NUM_THREADS=<n>. On an n-core ",
                     "machine, -t n.")
    end
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
