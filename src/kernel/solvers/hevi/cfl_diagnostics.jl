#=============================================================================
 cfl_diagnostics.jl -- direction-wise stability-limit report.

 WHY THIS EXISTS
 ---------------
 "The run is too slow, so make Δt bigger" is only actionable once you know
 which physical process is actually setting the explicit stability limit.
 On a stratified LES there are four candidates and they scale differently
 with the mesh:

     acoustic     Δ / (|u| + c)      c ≈ 340 m/s          ~ Δ
     advective    Δ / |u|            |u| ≈ 10 m/s         ~ Δ
     viscous      Δ² / (2 ν_sgs)     ν_t from the SGS     ~ Δ²
     buoyancy     1 / N              N ≈ 0.01 1/s         mesh-independent

 The viscous one is the trap: it goes as Δ², so a mesh refined vertically to
 resolve a surface layer can be limited by the SGS diffusion rather than by
 sound, and in that case a horizontally-explicit/vertically-implicit (HEVI)
 acoustic split buys nothing -- the cheap fix is to make the *vertical
 diffusion* implicit instead (see :hevi_implicit_diffusion in hevi/README.md).

 Rather than the usual Δ/N² rule of thumb for the LGL node clustering, this
 walks the actual element node graph and measures the true node-to-node
 distance along each reference direction. For nop=4 the Δ/N² estimate is
 pessimistic by ~2.8x (LGL min gap is 0.1727·h, not h/16), which is enough
 to change which term looks binding.

 USAGE
 -----
 From a deck:                 :lcfl_report => true
 Fires once, after the first RHS evaluation of the run (the SGS viscosity
 cache has to be populated before ν_t means anything), and again at every
 diagnostics output time when :lcfl_report_every is also true.

 Standalone, from a restart or the REPL:
     cfl_report(params, u, t)
=============================================================================#

"""
    NodeSpacings

Per-direction physical node spacing on the local mesh, in metres.

`hx`, `hy`, `hz` are indexed `[iel, i, j, k]` and hold the distance from node
`(i,j,k)` to its nearest neighbour along the ξ, η and ζ reference directions
respectively. They are *physical* distances (the full 3-vector norm), so they
stay meaningful on a warped or stretched mesh where the reference directions
are not aligned with x/y/z.
"""
struct NodeSpacings{A4 <: AbstractArray{<:AbstractFloat,4}}
    hx::A4
    hy::A4
    hz::A4
end

"""
    compute_node_spacings(mesh, coords) -> NodeSpacings

Measure the node-to-node distance along each reference direction, at every
node of every local element.

For an interior node this is the *minimum* of the distances to the two
neighbours along that direction -- the stability limit is set by the smallest
gap, and on an LGL grid the two gaps straddling a node differ by up to 2x.
Edge nodes (i == 1 or i == ngl) have only one neighbour along that direction,
so that single distance is used.
"""
function compute_node_spacings(mesh, coords)

    ngl   = Int(mesh.ngl)
    nelem = Int(mesh.nelem)
    T     = eltype(coords)

    hx = zeros(T, nelem, ngl, ngl, ngl)
    hy = zeros(T, nelem, ngl, ngl, ngl)
    hz = zeros(T, nelem, ngl, ngl, ngl)

    @inline function _dist(ipa, ipb)
        dx = coords[1, ipa] - coords[1, ipb]
        dy = coords[2, ipa] - coords[2, ipb]
        dz = coords[3, ipa] - coords[3, ipb]
        return sqrt(dx*dx + dy*dy + dz*dz)
    end

    conn = mesh.connijk
    @inbounds for iel = 1:nelem, k = 1:ngl, j = 1:ngl, i = 1:ngl
        ip = conn[iel, i, j, k]

        # ξ direction
        d = typemax(T)
        i > 1   && (d = min(d, _dist(ip, conn[iel, i-1, j, k])))
        i < ngl && (d = min(d, _dist(ip, conn[iel, i+1, j, k])))
        hx[iel, i, j, k] = d

        # η direction
        d = typemax(T)
        j > 1   && (d = min(d, _dist(ip, conn[iel, i, j-1, k])))
        j < ngl && (d = min(d, _dist(ip, conn[iel, i, j+1, k])))
        hy[iel, i, j, k] = d

        # ζ direction
        d = typemax(T)
        k > 1   && (d = min(d, _dist(ip, conn[iel, i, j, k-1])))
        k < ngl && (d = min(d, _dist(ip, conn[iel, i, j, k+1])))
        hz[iel, i, j, k] = d
    end

    return NodeSpacings(hx, hy, hz)
end

"""
    cfl_limits(params, u) -> NamedTuple

Compute the per-direction explicit stability limits over the *global* mesh
(MPI-reduced), in seconds. Every entry is the largest Δt for which the
corresponding CFL-style number equals one, so the number you can actually run
is this divided by the scheme's stability radius.

Returned fields (all `Float64`, all MPI `MIN`-reduced):

  `dt_acoustic_x/y/z`  h / (|u| + c)      sound in each reference direction
  `dt_advective_x/y/z` h / |u|            transport only
  `dt_viscous_x/y/z`   h² / (2 ν_eff)     SGS diffusion, parabolic limit
  `dt_buoyancy`        1 / N              Brunt-Vaisala, if N² > 0
  `hmin_x/y/z`         smallest node spacing in each direction, metres
  `cmax`, `umax`, `wmax`, `numax`  extrema of the fields that drive the above

`ν_eff` is the largest kinematic diffusivity actually applied by the SGS
closure, i.e. it includes the per-equation `:μ` mask and the turbulent
Prandtl number, because a scalar equation with `μ[5] = 2` is diffused twice
as hard as momentum and it is the largest one that binds Δt.
"""
function cfl_limits(params, u)

    mesh   = params.mesh
    ngl    = Int(mesh.ngl)
    nelem  = Int(mesh.nelem)
    neqs   = Int(params.neqs)
    coords = mesh.coords
    PhysConst = PhysicalConst{Float64}()
    γ  = PhysConst.γ
    g  = PhysConst.g

    # uaux must hold the current state: every rhs! call refreshes it from u,
    # but a caller may hand us a `u` that has never been through rhs!.
    u2uaux!(@view(params.uaux[:,:]), u, neqs, mesh.npoin)
    uaux = params.uaux

    sp = compute_node_spacings(mesh, coords)

    # Effective kinematic diffusivity per node.
    #
    # :μ means two different things depending on the closure, and reading it
    # the wrong way makes the viscous row of this table wrong in the direction
    # that matters -- it reports "inf", i.e. "diffusion never limits Δt".
    #
    #   AV()                :μ entries ARE diffusivities, in m²/s. There is no
    #                       μ_turb array at all (params.sgs is nothing).
    #   SMAG/VREM/DSGS      :μ entries are per-equation MULTIPLIERS on the
    #                       eddy viscosity the closure computed; the physical
    #                       quantity is μ_turb, which is DYNAMIC (ρ ν_t) -- see
    #                       the SGS_diffusion comment in SGS.jl -- so it is
    #                       divided by ρ below.
    have_sgs = params.inputs[:lvisc] == true &&
               hasproperty(params, :sgs) &&
               params.sgs !== nothing &&
               hasproperty(params.sgs, :μ_turb) &&
               length(params.sgs.μ_turb) == mesh.npoin
    μmask = params.visc_coeff
    is_av = params.inputs[:lvisc] == true && hasproperty(params, :VT) && params.VT isa AV
    ν_const = is_av ? maximum(Float64, μmask) : 0.0
    # momentum slot 2, scalar slot neqs: the two the closure treats differently
    mask_mom = length(μmask) >= 2 ? Float64(μmask[2]) : 1.0
    mask_sca = length(μmask) >= neqs ? Float64(μmask[neqs]) : 1.0
    Pr_t     = have_sgs ? Float64(params.sgs.Pr_t) : PhysConst.Pr_t
    μ_mol    = PhysConst.μ_mol

    BIG = 1.0e30
    dt_ax = BIG; dt_ay = BIG; dt_az = BIG
    dt_ux = BIG; dt_uy = BIG; dt_uz = BIG
    dt_vx = BIG; dt_vy = BIG; dt_vz = BIG
    hmin_x = BIG; hmin_y = BIG; hmin_z = BIG
    cmax = 0.0; umax = 0.0; wmax = 0.0; numax = 0.0
    N2max = 0.0

    conn = mesh.connijk
    @inbounds for iel = 1:nelem, k = 1:ngl, j = 1:ngl, i = 1:ngl
        ip = conn[iel, i, j, k]

        ρ  = uaux[ip, 1]
        ρ <= 0.0 && continue
        uu = uaux[ip, 2] / ρ
        vv = uaux[ip, 3] / ρ
        ww = uaux[ip, 4] / ρ
        θ  = uaux[ip, 5] / ρ
        p  = perfectGasLaw_ρθtoP(PhysConst, ρ, θ)
        c  = sqrt(max(γ * p / ρ, 0.0))

        ν = ν_const
        if have_sgs
            μt = params.sgs.μ_turb[ip]
            ν  = max(ν, max((μ_mol + μt) * mask_mom, (μt / Pr_t) * mask_sca) / ρ)
            if hasproperty(params.sgs, :N2)
                N2max = max(N2max, params.sgs.N2[ip])
            end
        end

        hx = sp.hx[iel, i, j, k]
        hy = sp.hy[iel, i, j, k]
        hz = sp.hz[iel, i, j, k]

        au = abs(uu); av = abs(vv); aw = abs(ww)
        cmax  = max(cmax, c)
        umax  = max(umax, max(au, av))
        wmax  = max(wmax, aw)
        numax = max(numax, ν)
        hmin_x = min(hmin_x, hx); hmin_y = min(hmin_y, hy); hmin_z = min(hmin_z, hz)

        dt_ax = min(dt_ax, hx / (au + c));  dt_ay = min(dt_ay, hy / (av + c))
        dt_az = min(dt_az, hz / (aw + c))
        au > 0 && (dt_ux = min(dt_ux, hx / au))
        av > 0 && (dt_uy = min(dt_uy, hy / av))
        aw > 0 && (dt_uz = min(dt_uz, hz / aw))
        if ν > 0
            dt_vx = min(dt_vx, hx*hx / (2ν)); dt_vy = min(dt_vy, hy*hy / (2ν))
            dt_vz = min(dt_vz, hz*hz / (2ν))
        end
    end

    comm = params.SD == NSD_1D() ? get_mpi_comm() : mesh.parts.comm
    _rmin(x) = MPI.Allreduce(x, MPI.MIN, comm)
    _rmax(x) = MPI.Allreduce(x, MPI.MAX, comm)

    dt_buoy = N2max > 0 ? 1.0 / sqrt(_rmax(N2max)) : BIG

    return (dt_acoustic_x  = _rmin(dt_ax), dt_acoustic_y  = _rmin(dt_ay), dt_acoustic_z  = _rmin(dt_az),
            dt_advective_x = _rmin(dt_ux), dt_advective_y = _rmin(dt_uy), dt_advective_z = _rmin(dt_uz),
            dt_viscous_x   = _rmin(dt_vx), dt_viscous_y   = _rmin(dt_vy), dt_viscous_z   = _rmin(dt_vz),
            dt_buoyancy    = dt_buoy,
            hmin_x = _rmin(hmin_x), hmin_y = _rmin(hmin_y), hmin_z = _rmin(hmin_z),
            cmax = _rmax(cmax), umax = _rmax(umax), wmax = _rmax(wmax), numax = _rmax(numax))
end

"""
    cfl_report(params, u, t; io = stdout, dt = nothing)

Print the direction-wise stability table on rank 0. `dt` defaults to the
deck's `:Δt` and is used to express each limit as the CFL number the run is
actually sitting at, which is the number you compare against 1.

The last block is the point of the whole exercise: it states what a HEVI
split would and would not remove, computed from these numbers rather than
from a rule of thumb.
"""
function cfl_report(params, u, t; io = stdout, dt = nothing)

    comm = params.SD == NSD_1D() ? get_mpi_comm() : params.mesh.parts.comm
    rank = MPI.Comm_rank(comm)
    L    = cfl_limits(params, u)
    rank == 0 || return L

    Δt = dt === nothing ? Float64(params.inputs[:Δt]) : Float64(dt)
    BIG = 1.0e29
    fmt(x) = x > BIG ? "      inf" : @sprintf("%9.4g", x)
    cfl(x) = x > BIG ? "     --  " : @sprintf("%9.4g", Δt / x)

    println(io)
    println(io, " ┌── CFL / stability report ───────────────────── t = ", @sprintf("%.4f", t), " s")
    println(io, " │  Δt in use = ", @sprintf("%.6g", Δt), " s")
    println(io, " │")
    println(io, " │  node spacing [m]     ξ:", fmt(L.hmin_x), "   η:", fmt(L.hmin_y), "   ζ:", fmt(L.hmin_z))
    @printf(io, " │  anisotropy  min(hx,hy)/hz = %.3f\n", min(L.hmin_x, L.hmin_y) / L.hmin_z)
    @printf(io, " │  extrema     c = %.1f m/s   |u_h| = %.2f m/s   |w| = %.3f m/s   nu_sgs = %.3g m^2/s\n",
            L.cmax, L.umax, L.wmax, L.numax)
    println(io, " │")
    println(io, " │  process            limit dt [s]                 CFL at this dt")
    println(io, " │                   xi        eta       zeta        xi        eta       zeta")
    println(io, " │  acoustic  ", fmt(L.dt_acoustic_x),  fmt(L.dt_acoustic_y),  fmt(L.dt_acoustic_z),
                                  cfl(L.dt_acoustic_x),  cfl(L.dt_acoustic_y),  cfl(L.dt_acoustic_z))
    println(io, " │  advective ", fmt(L.dt_advective_x), fmt(L.dt_advective_y), fmt(L.dt_advective_z),
                                  cfl(L.dt_advective_x), cfl(L.dt_advective_y), cfl(L.dt_advective_z))
    println(io, " │  viscous   ", fmt(L.dt_viscous_x),   fmt(L.dt_viscous_y),   fmt(L.dt_viscous_z),
                                  cfl(L.dt_viscous_x),   cfl(L.dt_viscous_y),   cfl(L.dt_viscous_z))
    println(io, " │  buoyancy  ", fmt(L.dt_buoyancy))
    println(io, " │")

    #-------------------------------------------------------------------------
    # What a HEVI split would buy, from these numbers.
    #
    # An explicit scheme is limited by the SUM of the per-direction rates
    # (the eigenvalues add along the diagonal of the tensor-product operator),
    # not by the smallest of them, so the comparison has to be made on rates.
    # HEVI removes exactly one term from that sum: the vertical acoustic one.
    # Everything else -- vertical advection, vertical diffusion, both
    # horizontal acoustics -- stays explicit.
    #-------------------------------------------------------------------------
    _rate(x) = x > BIG ? 0.0 : 1.0 / x
    rate_full = _rate(L.dt_acoustic_x) + _rate(L.dt_acoustic_y) + _rate(L.dt_acoustic_z) +
                _rate(L.dt_viscous_x)  + _rate(L.dt_viscous_y)  + _rate(L.dt_viscous_z)
    rate_hevi = _rate(L.dt_acoustic_x) + _rate(L.dt_acoustic_y) + _rate(L.dt_advective_z) +
                _rate(L.dt_viscous_x)  + _rate(L.dt_viscous_y)  + _rate(L.dt_viscous_z)
    rate_hevi_diff = _rate(L.dt_acoustic_x) + _rate(L.dt_acoustic_y) + _rate(L.dt_advective_z) +
                     _rate(L.dt_viscous_x)  + _rate(L.dt_viscous_y)
    # substepping the acoustics on top of HEVI: outer step sees no sound at all
    rate_split = _rate(L.dt_advective_x) + _rate(L.dt_advective_y) + _rate(L.dt_advective_z) +
                 _rate(L.dt_viscous_x)   + _rate(L.dt_viscous_y)   + _rate(L.dt_viscous_z)
    # ... and with the vertical diffusion implicit as well (:implicit_vdiff),
    # which is the only combination that leaves NOTHING vertical in the
    # explicit budget. On a mesh refined in z this is usually the largest
    # number in the table, and on a spun-up boundary layer it is often the only
    # one that differs from the row above it.
    rate_split_diff = _rate(L.dt_advective_x) + _rate(L.dt_advective_y) +
                      _rate(L.dt_advective_z) +
                      _rate(L.dt_viscous_x)   + _rate(L.dt_viscous_y)

    gain(r) = r > 0 ? rate_full / r : Inf
    println(io, " │  Δt gain available, per scheme (rate-summed, before per-stage cost):")
    @printf(io, " │    HEVI  (vertical acoustics implicit)            %6.2fx\n", gain(rate_hevi))
    @printf(io, " │    HEVI + implicit vertical diffusion            %6.2fx\n", gain(rate_hevi_diff))
    @printf(io, " │    HEVI + acoustic substepping (outer step)      %6.2fx\n", gain(rate_split))
    # IMEX3D removes the same three acoustic terms the substepping outer step
    # does, so the Δt gain is identical. What differs is the price: substepping
    # pays it in inner steps at the horizontal acoustic limit, IMEX3D in Krylov
    # iterations per stage. Listed separately because they are separate schemes
    # and the reader is choosing between them.
    @printf(io, " │    IMEX3D (all acoustics implicit)               %6.2fx\n", gain(rate_split))
    @printf(io, " │    IMEX3D + implicit vertical diffusion          %6.2fx\n", gain(rate_split_diff))
    println(io, " │")
    binder = argmax([_rate(L.dt_acoustic_x) + _rate(L.dt_acoustic_y),
                     _rate(L.dt_acoustic_z),
                     _rate(L.dt_viscous_x) + _rate(L.dt_viscous_y) + _rate(L.dt_viscous_z),
                     _rate(L.dt_advective_x) + _rate(L.dt_advective_y) + _rate(L.dt_advective_z)])
    names = ("horizontal acoustics -- HEVI cannot help; substepping or IMEX3D removes it",
             "vertical acoustics   -- this is what HEVI is for",
             "SGS diffusion        -- :implicit_vdiff => true; the acoustics are not what binds",
             "advection            -- already at the floor; no time-integration trick helps")
    println(io, " │  dominant term: ", names[binder])
    println(io, " └", "─"^70)
    flush(io)

    return L
end
