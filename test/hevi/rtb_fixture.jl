#=============================================================================
 rtb_fixture.jl -- a rising thermal bubble on the mock spectral-element mesh.

 The classic dry nonhydrostatic benchmark (Robert 1993; Giraldo & Restelli,
 J. Comput. Phys. 227, 2008): a warm, smooth perturbation in a neutrally
 stratified atmosphere at rest, which accelerates upward under its own
 buoyancy and deforms into a mushroom.

 WHY THIS IS THE RIGHT CASE FOR HEVI
 -----------------------------------
 It exercises precisely the three things the split touches, and nothing it
 does not:

   * vertical acoustics -- the terms HEVI moves to the implicit side;
   * buoyancy           -- the -g dρ coupling in the implicit operator, and
                           the discrete hydrostatic balance it must not break;
   * a base state at rest -- so any spurious motion is visible. A split that
     quietly disturbed the balance would show up as the whole domain drifting,
     not as a subtle error in the bubble.

 The RHS here is the full nonlinear compressible Euler system in (ρ, ρu, ρv,
 ρw, ρθ) with gravity, discretised exactly as Jexpresso discretises it:
 strong form, collocation differentiation, -ωJ(div - S) into the element RHS,
 then DSS, the MPI assembly and the inverse mass matrix. It is a fixture, not
 a second solver -- it exists so the HEVI machinery can be driven by real
 buoyancy-driven nonlinear dynamics without the full dependency stack.

 Walls are free-slip and are imposed the way Jexpresso imposes them: by
 projecting the normal momentum out of the STATE at the start of the RHS, not
 by zeroing du. That matters for the split. With the projection on the state,
 f_imp is evaluated on the same projected state the flux divergence saw, so
 f_exp = f_full - f_imp still holds exactly at the wall. Zeroing du instead
 would put the condition in f_exp only, and the wall momentum would drift
 under any tableau whose explicit and implicit weights differ.
=============================================================================#

"""
    rtb_initial_state(params; θc, r0, xc, yc, zc) -> u0

Neutrally stratified hydrostatic base state plus a smooth warm bubble.

The perturbation is the cosine bump

    Δθ = θc/2 · (1 + cos(π r / r0))    for r < r0

rather than the linear taper `θc(1 - r/r0)`. Both appear in the literature;
the cosine one is C¹ at the bubble edge, and on a spectral element mesh a
kink there rings for the rest of the run. This is a test of the time
integrator, so the spatial solution should not be fighting Gibbs oscillations
at the same time.

`r` is measured in the x–z plane only, so the bubble is a cylinder along y and
the solution stays y-invariant — which makes "did anything drift in y?" a
usable diagnostic.
"""
function rtb_initial_state(params; θc = 0.5, r0 = 250.0,
                           xc = 500.0, zc = 350.0, θ0 = 300.0)

    PC    = PhysicalConst{Float64}()
    mesh  = params.mesh
    npoin = mesh.npoin
    qe    = params.qp.qe
    u0    = zeros(npoin * 5)

    for ip = 1:npoin
        x, z = mesh.coords[1, ip], mesh.coords[3, ip]
        r  = sqrt((x - xc)^2 + (z - zc)^2)
        Δθ = r < r0 ? 0.5 * θc * (1.0 + cospi(r / r0)) : 0.0
        θ  = θ0 + Δθ

        # Hydrostatic pressure of the BASE state, with the perturbation carried
        # by θ alone: the bubble is buoyant, not over-pressured.
        pbar = PC.pref * (1.0 - PC.g * z / (PC.cp * θ0))^(PC.cp / PC.Rair)
        ρ    = (1.0 / θ) * (pbar / PC.C0)^(1.0 / PC.γ)

        u0[0*npoin + ip] = ρ
        u0[1*npoin + ip] = 0.0
        u0[2*npoin + ip] = 0.0
        u0[3*npoin + ip] = 0.0
        u0[4*npoin + ip] = ρ * θ
    end
    return u0
end

"""
    RTBWalls(mesh) -> NamedTuple of boundary node lists

Nodes on each pair of faces of the box, for the free-slip projection.
"""
function rtb_walls(mesh)
    tol = 1.0e-6 * max(mesh.xmax - mesh.xmin, mesh.zmax - mesh.zmin)
    xw = Int[]; yw = Int[]; zw = Int[]
    for ip = 1:mesh.npoin
        x, y, z = mesh.coords[1, ip], mesh.coords[2, ip], mesh.coords[3, ip]
        (abs(x - mesh.xmin) < tol || abs(x - mesh.xmax) < tol) && push!(xw, ip)
        (abs(y - mesh.ymin) < tol || abs(y - mesh.ymax) < tol) && push!(yw, ip)
        (abs(z - mesh.zmin) < tol || abs(z - mesh.zmax) < tol) && push!(zw, ip)
    end
    return (x = xw, y = yw, z = zw)
end

"""
    RTBCache

Scratch for the Euler RHS: element-local fluxes and source, the element RHS,
the assembled nodal RHS and its MPI assembler.
"""
struct RTBCache
    F::Array{Float64,4}; G::Array{Float64,4}; H::Array{Float64,4}; S::Array{Float64,4}
    rhs_el::Array{Float64,5}
    RHS::Matrix{Float64}
    vaux::Vector{Float64}
    dss::MockAssembler
    walls::NamedTuple{(:x, :y, :z), Tuple{Vector{Int}, Vector{Int}, Vector{Int}}}
end

function RTBCache(params)
    ngl   = Int(params.mesh.ngl)
    nelem = Int(params.mesh.nelem)
    npoin = Int(params.mesh.npoin)
    RHS   = zeros(npoin, 5)
    RTBCache(zeros(ngl, ngl, ngl, 5), zeros(ngl, ngl, ngl, 5),
             zeros(ngl, ngl, ngl, 5), zeros(ngl, ngl, ngl, 5),
             zeros(nelem, ngl, ngl, ngl, 5), RHS, zeros(npoin),
             setup_assembler(NSD_3D(), RHS, params.mesh.ip2gip, params.mesh.gip2owner),
             rtb_walls(params.mesh))
end

"""
    rtb_rhs!(du, u, p, t)

Full nonlinear compressible Euler with gravity, the "explicit RHS" the ARK
integrator is handed. `p` must carry `p.rtb::RTBCache` and, for a HEVI run,
`p.hevi`.

Mutates `u`: the free-slip projection is applied to the state, as described at
the top of this file.
"""
function rtb_rhs!(du, u, p, t)

    params = p
    mesh   = params.mesh
    cache  = p.rtb
    ngl    = Int(mesh.ngl)
    nelem  = Int(mesh.nelem)
    npoin  = Int(mesh.npoin)
    conn   = mesh.connijk
    dψ     = params.basis.dψ
    ω      = params.ω
    met    = params.metrics
    PC     = PhysicalConst{Float64}()
    γ, C0, g = PC.γ, PC.C0, PC.g

    # free-slip walls, on the state
    @inbounds for ip in cache.walls.x; u[1*npoin + ip] = 0.0; end
    @inbounds for ip in cache.walls.y; u[2*npoin + ip] = 0.0; end
    @inbounds for ip in cache.walls.z; u[3*npoin + ip] = 0.0; end

    F, G, H, S = cache.F, cache.G, cache.H, cache.S
    rhs_el = cache.rhs_el
    fill!(rhs_el, 0.0)

    @inbounds for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = conn[iel, i, j, k]
            ρ  = u[0*npoin + ip]; ρu = u[1*npoin + ip]; ρv = u[2*npoin + ip]
            ρw = u[3*npoin + ip]; ρθ = u[4*npoin + ip]
            uu = ρu/ρ; vv = ρv/ρ; ww = ρw/ρ
            # Clamped so a diverging solution yields Inf/NaN instead of
            # throwing. (ρθ)^γ on a negative ρθ raises DomainError, and an
            # exception raised inside the RHS is a DEADLOCK under MPI: the
            # throwing rank unwinds out of assemble_mpi! while the others sit
            # in it. Divergence has to stay detectable by isfinite, not by
            # catching. The clamp is exact wherever the solution is valid.
            P  = C0 * max(ρθ, 0.0)^γ

            F[i,j,k,1] = ρu;      F[i,j,k,2] = ρu*uu + P; F[i,j,k,3] = ρu*vv
            F[i,j,k,4] = ρu*ww;   F[i,j,k,5] = ρθ*uu
            G[i,j,k,1] = ρv;      G[i,j,k,2] = ρv*uu;     G[i,j,k,3] = ρv*vv + P
            G[i,j,k,4] = ρv*ww;   G[i,j,k,5] = ρθ*vv
            H[i,j,k,1] = ρw;      H[i,j,k,2] = ρw*uu;     H[i,j,k,3] = ρw*vv
            H[i,j,k,4] = ρw*ww + P; H[i,j,k,5] = ρθ*ww

            S[i,j,k,1] = 0.0; S[i,j,k,2] = 0.0; S[i,j,k,3] = 0.0
            S[i,j,k,4] = -ρ*g
            S[i,j,k,5] = 0.0
        end

        for ieq = 1:5, k = 1:ngl, j = 1:ngl, i = 1:ngl
            dFdξ = 0.0; dGdη = 0.0; dHdζ = 0.0
            for m = 1:ngl
                dFdξ += dψ[m,i] * F[m,j,k,ieq]
                dGdη += dψ[m,j] * G[i,m,k,ieq]
                dHdζ += dψ[m,k] * H[i,j,m,ieq]
            end
            div  = dFdξ*met.dξdx[iel,i,j,k] + dGdη*met.dηdy[iel,i,j,k] +
                   dHdζ*met.dζdz[iel,i,j,k]
            ωJac = ω[i]*ω[j]*ω[k]*met.Je[iel,i,j,k]
            rhs_el[iel,i,j,k,ieq] -= ωJac * (div - S[i,j,k,ieq])
        end
    end

    RHS = cache.RHS
    fill!(RHS, 0.0)
    DSS_rhs!(RHS, rhs_el, conn, nelem, ngl, 5, NSD_3D(), ContGal())
    assemble_mpi!(RHS, cache.dss)
    for ieq = 1:5
        divide_by_mass_matrix!(@view(RHS[:, ieq]), cache.vaux, params.Minv, 5, npoin, ContGal())
    end

    @inbounds for ieq = 1:5
        off = (ieq-1)*npoin
        for ip = 1:npoin
            du[off + ip] = RHS[ip, ieq]
        end
    end
    return du
end

"""
    rtb_diagnostics(u, params) -> NamedTuple

`θmax`  largest θ perturbation above the 300 K base, K
`zbar`  θ'-weighted mean height of the perturbation, m -- the bubble's centre
        of buoyancy, and the thing that must go UP
`wmax`  largest |w|, m/s
`vmax`  largest |v|, m/s -- the bubble is a cylinder along y, so this stays at
        round-off unless something has gone wrong
"""
function rtb_diagnostics(u, params; θ0 = 300.0)
    mesh  = params.mesh
    npoin = mesh.npoin
    rank  = MPI.Comm_rank(mesh.parts.comm)
    θmax = -Inf; wmax = 0.0; vmax = 0.0
    num = 0.0; den = 0.0
    for ip = 1:npoin
        ρ  = u[0*npoin + ip]
        θp = u[4*npoin + ip]/ρ - θ0
        w  = u[3*npoin + ip]/ρ
        v  = u[2*npoin + ip]/ρ
        θmax = max(θmax, θp); wmax = max(wmax, abs(w)); vmax = max(vmax, abs(v))
        # The weighted mean is a SUM, so it has to run over owned nodes only:
        # a node on a rank interface exists on both ranks and would otherwise
        # be counted twice, making zbar depend on the rank count. The maxima
        # above are immune -- a duplicate cannot change a maximum.
        if θp > 0.01 && mesh.gip2owner[ip] == rank
            num += θp * mesh.coords[3, ip]; den += θp
        end
    end
    comm = mesh.parts.comm
    if MPI.Comm_size(comm) > 1
        θmax = MPI.Allreduce(θmax, MPI.MAX, comm)
        wmax = MPI.Allreduce(wmax, MPI.MAX, comm)
        vmax = MPI.Allreduce(vmax, MPI.MAX, comm)
        num  = MPI.Allreduce(num,  MPI.SUM, comm)
        den  = MPI.Allreduce(den,  MPI.SUM, comm)
    end
    return (θmax = θmax, zbar = den > 0 ? num/den : NaN, wmax = wmax, vmax = vmax)
end
