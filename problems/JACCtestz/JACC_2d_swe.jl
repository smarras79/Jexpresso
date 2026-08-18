#---------------------------------------------------------------------------------
#
# JACC_2d_swe.jl — a STANDALONE exercise of the portable 2D shallow-water RHS in
#                  src/kernel/operators/rhs_jacc.jl.
#
# The 2D companion of problems/JACCtestz/JACC_1d_gpu.jl: it needs JACC and nothing
# else — no MPI, no Gridap, no gmsh, no Jexpresso module — so it runs in seconds
# and is the first thing to try when a JACC backend is being brought up on a new
# machine.
#
#   julia --project -t 8 problems/JACCtestz/JACC_2d_swe.jl
#   julia --project -t 8 problems/JACCtestz/JACC_2d_swe.jl 12 8 4     # nelx nely nop
#
# For a GPU run, set the backend ONCE (it is a Preferences value, so it needs a
# restart) and then run exactly the same command:
#
#   julia --project -e 'using JACC; JACC.set_backend("cuda")'   # or amdgpu/oneapi/metal
#   julia --project problems/JACCtestz/JACC_2d_swe.jl
#
# WHAT IT CHECKS
# --------------
#   1. well-balancedness   the lake at rest is an exact steady state of the
#                          discrete operator: RHS(qe) must be zero to round-off.
#                          This is a real check on the flux/source split of
#                          problems/ShallowWater/SoliWaveIsland/user_flux.jl, not
#                          just on the plumbing.
#   2. correctness         the JACC RHS against an independent serial reference
#                          written straight from the weak form, on the solitary
#                          wave initial state, with and without viscosity.
#   3. determinism         two runs agree BIT FOR BIT. They do because the direct
#                          stiffness summation is a gather over a fixed CSR map,
#                          not an atomic scatter — see the header of rhs_jacc.jl.
#   4. throughput          RHS evaluations per second on the active backend.
#
# The mesh is built here, curvilinearly (the interior nodes are warped, the four
# walls stay straight), so the metric terms are non-trivial and dξdy / dηdx are
# genuinely non-zero — a Cartesian test would let an index slip through.
#
#---------------------------------------------------------------------------------

using JACC
using Printf
using LinearAlgebra

#---------------------------------------------------------------------------------
# The four device callbacks of the SoliWaveIsland deck.
#
# Copied verbatim from problems/ShallowWater/SoliWaveIsland/{user_flux,user_source,
# user_primitives,user_bc}.jl. They are copied and not `include`d because those
# files ALSO define the host `user_flux!(…, ::St_mesh, ::CL, ::TOTAL)` methods,
# whose signatures name Jexpresso types that do not exist in this standalone
# script. Keep the two in sync; the physics lives in the deck, this is a mirror.
#---------------------------------------------------------------------------------
const _G_SWE      = 9.81
const _H_WET_SWE  = 1.0e-3

const _XC_CONE_SWE   =  12.5
const _YC_CONE_SWE   =   0.0
const _RC_CONE_SWE   =   3.6
const _HC_CONE_SWE   =   0.93
const _SIGMA_DRY_SWE =  25.0

function user_flux_gpu(q, qe, PhysConst, lpert)
    T  = eltype(q)
    Hc = max(q[1], T(0.0))
    He = qe[1]
    ε  = T(_H_WET_SWE)
    H4 = max(Hc, ε)
    de = sqrt(T(2.0)) * Hc / sqrt(Hc^4 + H4^4)
    u  = de * q[2]
    v  = de * q[3]
    p  = T(0.5) * T(_G_SWE) * (Hc * Hc - He * He)

    return T(Hc*u), T(Hc*u*u + p), T(Hc*u*v),
           T(Hc*v), T(Hc*v*u),     T(Hc*v*v + p)
end

function user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)
    T  = eltype(q)
    H  = q[1]
    dH = max(H, T(0.0)) - qe[1]

    dx = x - T(_XC_CONE_SWE)
    dy = y - T(_YC_CONE_SWE)
    r  = sqrt(dx*dx + dy*dy)

    dHbdx = T(0.0)
    dHbdy = T(0.0)
    if r < T(_RC_CONE_SWE) && r > T(1.0e-12)
        slope = -T(_HC_CONE_SWE) / (T(_RC_CONE_SWE) * r)
        dHbdx = slope * dx
        dHbdy = slope * dy
    end

    g = T(_G_SWE)
    σ = H < T(_H_WET_SWE) ? T(_SIGMA_DRY_SWE) : T(0.0)
    return T(0.0), T(-g * dH * dHbdx - σ * q[2]), T(-g * dH * dHbdy - σ * q[3])
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    return T(u[1] - qe[1]), T(u[2]), T(u[3])
end

function user_bc_dirichlet_gpu(q, qe, coords, t, nx, ny, qbdy, lpert)
    x, y = coords[1], coords[2]
    T  = eltype(q)
    Hu = q[2]
    Hv = q[3]
    qn = nx * Hu + ny * Hv
    u  = Hu - qn * nx
    v  = Hv - qn * ny
    return T(qbdy[1]), T(u), T(v)
end

# The kernels under test. Only the kernel + cache definitions are used; the
# `params`-driven entry points at the bottom of that file are Jexpresso's.
include(joinpath(@__DIR__, "..", "..", "src", "kernel", "operators", "rhs_jacc.jl"))


#=================================================================================
                     A SMALL CURVILINEAR SEM MESH, BUILT HERE
=================================================================================#

#--- Legendre-Gauss-Lobatto nodes and weights on [-1,1] --------------------------
function lgl_nodes_weights(N::Int)
    ngl = N + 1
    ξ = zeros(Float64, ngl)
    w = zeros(Float64, ngl)

    # Chebyshev-Gauss-Lobatto start, then Newton on (1-x²)P'_N(x)
    for i = 0:N
        x = -cos(pi*i/N)
        for _ = 1:100
            # Legendre P_N and P_{N-1} by recurrence
            p0, p1 = 1.0, x
            for k = 2:N
                p0, p1 = p1, ((2k-1)*x*p1 - (k-1)*p0)/k
            end
            PN, PNm1 = p1, p0
            dPN  = N*(x*PN - PNm1)/(x*x - 1 + eps())
            d2PN = (2x*dPN - N*(N+1)*PN)/(1 - x*x + eps())
            abs(dPN) < 1e-300 && break
            dx = -dPN/d2PN
            x += dx
            abs(dx) < 1e-15 && break
        end
        ξ[i+1] = x
    end
    ξ[1], ξ[end] = -1.0, 1.0
    ξ[abs.(ξ) .< 1e-15] .= 0.0

    for i = 1:ngl
        x = ξ[i]
        p0, p1 = 1.0, x
        for k = 2:N
            p0, p1 = p1, ((2k-1)*x*p1 - (k-1)*p0)/k
        end
        PN = N == 0 ? 1.0 : p1
        w[i] = 2.0/(N*(N+1)*PN*PN)
    end
    return ξ, w
end

#--- dψ[a,b] = ℓ_a'(ξ_b), the nodal derivative matrix on the LGL points -----------
function lgl_dpsi(ξ::Vector{Float64})
    ngl = length(ξ)
    N   = ngl - 1
    L   = similar(ξ)                       # P_N at the nodes
    for i = 1:ngl
        x = ξ[i]
        p0, p1 = 1.0, x
        for k = 2:N
            p0, p1 = p1, ((2k-1)*x*p1 - (k-1)*p0)/k
        end
        L[i] = p1
    end

    dψ = zeros(Float64, ngl, ngl)
    for a = 1:ngl, b = 1:ngl
        if a != b
            dψ[a,b] = (L[b]/L[a]) / (ξ[b] - ξ[a])
        end
    end
    dψ[1,1]     = -N*(N+1)/4
    dψ[ngl,ngl] =  N*(N+1)/4
    return dψ
end

#--- the physical map: a rectangle whose INTERIOR is warped, walls left straight --
@inline function _warp(s, t, Lx, Ly, A)
    b = A*sin(pi*s)*sin(pi*t)
    return Lx*s + b, Ly*t + b
end

struct TestMesh
    npoin::Int; nelem::Int; ngl::Int; nelx::Int; nely::Int
    x::Vector{Float64}; y::Vector{Float64}
    coords::Matrix{Float64}      # (2, npoin) — as St_mesh stores it
    connijk::Array{Int,3}
    dξdx::Array{Float64,3}; dξdy::Array{Float64,3}
    dηdx::Array{Float64,3}; dηdy::Array{Float64,3}
    Je::Array{Float64,3}
    dψ::Matrix{Float64}; ω::Vector{Float64}
    Minv::Vector{Float64}
    poin_in_bdy_edge::Matrix{Int}
    nx::Matrix{Float64}; ny::Matrix{Float64}
    nedges_bdy::Int
    xmin::Float64; xmax::Float64; ymin::Float64; ymax::Float64
end

function build_test_mesh(nelx::Int, nely::Int, nop::Int;
                         Lx = 25.0, Ly = 10.0, A = 0.35)

    ngl = nop + 1
    ξ, ω = lgl_nodes_weights(nop)
    dψ   = lgl_dpsi(ξ)

    npx = nelx*nop + 1
    npy = nely*nop + 1
    npoin = npx*npy
    nelem = nelx*nely

    # global node coordinates on the reference [0,1]² lattice, then warped
    x = zeros(npoin); y = zeros(npoin)
    sgrid = zeros(npx); tgrid = zeros(npy)
    for ex = 1:nelx, i = 1:ngl
        sgrid[(ex-1)*nop + i] = ((ex-1) + (ξ[i]+1)/2)/nelx
    end
    for ey = 1:nely, j = 1:ngl
        tgrid[(ey-1)*nop + j] = ((ey-1) + (ξ[j]+1)/2)/nely
    end
    for jy = 1:npy, ix = 1:npx
        ip = (jy-1)*npx + ix
        x[ip], y[ip] = _warp(sgrid[ix], tgrid[jy], Lx, Ly, A)
    end
    # shift y so the cone centre (12.5, 0) of the deck sits inside the basin
    y .-= Ly/2

    connijk = zeros(Int, nelem, ngl, ngl)
    for ey = 1:nely, ex = 1:nelx
        iel = (ey-1)*nelx + ex
        for j = 1:ngl, i = 1:ngl
            connijk[iel,i,j] = ((ey-1)*nop + j - 1)*npx + (ex-1)*nop + i
        end
    end

    # metric terms, exactly as metric_terms.jl does them: differentiate the nodal
    # coordinates with dψ, then invert the 2x2 Jacobian node by node.
    dξdx = zeros(nelem, ngl, ngl); dξdy = zeros(nelem, ngl, ngl)
    dηdx = zeros(nelem, ngl, ngl); dηdy = zeros(nelem, ngl, ngl)
    Je   = zeros(nelem, ngl, ngl)
    for iel = 1:nelem, j = 1:ngl, i = 1:ngl
        dxdξ = dxdη = dydξ = dydη = 0.0
        for k = 1:ngl
            dxdξ += dψ[k,i]*x[connijk[iel,k,j]]
            dxdη += dψ[k,j]*x[connijk[iel,i,k]]
            dydξ += dψ[k,i]*y[connijk[iel,k,j]]
            dydη += dψ[k,j]*y[connijk[iel,i,k]]
        end
        J = dxdξ*dydη - dxdη*dydξ
        J > 0 || error("non-positive Jacobian at element $iel node ($i,$j): J = $J")
        Je[iel,i,j]   = J
        dξdx[iel,i,j] =  dydη/J
        dξdy[iel,i,j] = -dxdη/J
        dηdx[iel,i,j] = -dydξ/J
        dηdy[iel,i,j] =  dxdξ/J
    end

    # diagonal CG mass matrix
    M = zeros(npoin)
    for iel = 1:nelem, j = 1:ngl, i = 1:ngl
        M[connijk[iel,i,j]] += ω[i]*ω[j]*Je[iel,i,j]
    end
    Minv = 1.0 ./ M

    # boundary edges of the rectangle, with their outward normals
    nedges_bdy = 2*nelx + 2*nely
    poin_in_bdy_edge = zeros(Int, nedges_bdy, ngl)
    nxb = zeros(nedges_bdy, ngl); nyb = zeros(nedges_bdy, ngl)
    e = 0
    for ex = 1:nelx                                  # south, n = (0,-1)
        e += 1
        for i = 1:ngl
            poin_in_bdy_edge[e,i] = connijk[ex, i, 1]
            nxb[e,i] = 0.0; nyb[e,i] = -1.0
        end
    end
    for ex = 1:nelx                                  # north, n = (0,+1)
        e += 1
        iel = (nely-1)*nelx + ex
        for i = 1:ngl
            poin_in_bdy_edge[e,i] = connijk[iel, i, ngl]
            nxb[e,i] = 0.0; nyb[e,i] = 1.0
        end
    end
    for ey = 1:nely                                  # west, n = (-1,0)
        e += 1
        iel = (ey-1)*nelx + 1
        for j = 1:ngl
            poin_in_bdy_edge[e,j] = connijk[iel, 1, j]
            nxb[e,j] = -1.0; nyb[e,j] = 0.0
        end
    end
    for ey = 1:nely                                  # east, n = (+1,0)
        e += 1
        iel = (ey-1)*nelx + nelx
        for j = 1:ngl
            poin_in_bdy_edge[e,j] = connijk[iel, ngl, j]
            nxb[e,j] = 1.0; nyb[e,j] = 0.0
        end
    end
    e == nedges_bdy || error("boundary edge count mismatch: $e vs $nedges_bdy")

    coords = permutedims(hcat(x, y))          # (2, npoin)

    return TestMesh(npoin, nelem, ngl, nelx, nely, x, y, coords, connijk,
                    dξdx, dξdy, dηdx, dηdy, Je, dψ, ω, Minv,
                    poin_in_bdy_edge, nxb, nyb, nedges_bdy,
                    minimum(x), maximum(x), minimum(y), maximum(y))
end


#=================================================================================
                        INITIAL STATE (the SoliWaveIsland deck)
=================================================================================#
function soliwave_state(m::TestMesh)
    neq = 3
    g   = _G_SWE
    A   = 0.064
    h0  = 0.32
    xcw = 2.5
    γ   = sqrt(3.0*A/(4.0*h0^3))
    H_dry = _H_WET_SWE

    q  = zeros(m.npoin, neq+1)
    qe = zeros(m.npoin, neq+1)
    for ip = 1:m.npoin
        xx = m.x[ip]; yy = m.y[ip]
        η  = A/cosh(γ*(xx - xcw))^2

        dx = xx - _XC_CONE_SWE; dy = yy - _YC_CONE_SWE
        r  = sqrt(dx*dx + dy*dy)
        Hb = r < _RC_CONE_SWE ? _HC_CONE_SWE*(1.0 - r/_RC_CONE_SWE) : 0.0

        H  = max(h0 + η - Hb, H_dry)
        ux = H > H_dry ? sqrt(g*(h0 + η))*η/(h0 + η) : 0.0

        q[ip,1] = H;  q[ip,2] = H*ux;  q[ip,3] = 0.0
        qe[ip,1] = max(h0 - Hb, H_dry)
    end
    return q, qe
end


#=================================================================================
    THE REFERENCE: the same weak form, written straight out, serial, no JACC
=================================================================================#
function reference_rhs(m::TestMesh, uaux::Matrix{Float64}, qe::Matrix{Float64},
                       t::Float64; neq = 3, lvisc = false, μ = Float64[],
                       lsource = true)

    ngl = m.ngl; nelem = m.nelem; npoin = m.npoin
    uaux = copy(uaux)

    # --- boundary values, edge by edge, with the CPU path's AlmostEqual gate ---
    # (Kopriva alg. 139, ε = 1e-6 — see _jacc_almost_equal in rhs_jacc.jl and
    #  AlmostEqual in src/kernel/infrastructure/Kopriva_functions.jl.)
    almost_equal(a, b) = (a == 0 || b == 0 || a <= 1e-6 || b <= 1e-6) ?
                         abs(a-b) <= 2e-6 :
                         (abs(a-b) <= 1e-6*abs(a) && abs(a-b) <= 1e-6*abs(b))
    for iedge = 1:m.nedges_bdy, k = 1:ngl
        ip = m.poin_in_bdy_edge[iedge,k]
        qb = fill(1234567.0, neq)
        vals = user_bc_dirichlet_gpu(view(uaux, ip, :), view(qe, ip, :),
                                     view(m.coords, :, ip), t,
                                     m.nx[iedge,k], m.ny[iedge,k], qb, false)
        for ieq = 1:neq
            v = vals[ieq]
            if !almost_equal(v, 1234567.0) && !almost_equal(v, uaux[ip,ieq])
                uaux[ip,ieq] = v
            end
        end
    end

    RHS = zeros(npoin, neq)
    F = zeros(ngl, ngl, 2neq); S = zeros(ngl, ngl, neq)
    U = zeros(ngl, ngl, neq)
    gξ = zeros(ngl, ngl, neq); gη = zeros(ngl, ngl, neq)

    for iel = 1:nelem
        for j = 1:ngl, i = 1:ngl
            ip = m.connijk[iel,i,j]
            F[i,j,:] .= user_flux_gpu(view(uaux, ip, 1:neq), view(qe, ip, :), nothing, false)
            if lsource
                S[i,j,:] .= user_source_gpu(view(uaux, ip, 1:neq), view(qe, ip, :),
                                            m.x[ip], m.y[ip], nothing,
                                            m.xmax, m.xmin, m.ymax, m.ymin, false)
            else
                S[i,j,:] .= 0.0
            end
            if lvisc
                U[i,j,:] .= user_primitives_gpu(view(uaux, ip, 1:neq), view(qe, ip, 1:neq), false)
            end
        end

        if lvisc
            for j = 1:ngl, i = 1:ngl
                ωJ = m.ω[i]*m.ω[j]*m.Je[iel,i,j]
                a11 = m.dξdx[iel,i,j]; a12 = m.dξdy[iel,i,j]
                a21 = m.dηdx[iel,i,j]; a22 = m.dηdy[iel,i,j]
                for ieq = 1:neq
                    dqdξ = 0.0; dqdη = 0.0
                    for k = 1:ngl
                        dqdξ += m.dψ[k,i]*U[k,j,ieq]
                        dqdη += m.dψ[k,j]*U[i,k,ieq]
                    end
                    qx = μ[ieq]*(dqdξ*a11 + dqdη*a21)
                    qy = μ[ieq]*(dqdξ*a12 + dqdη*a22)
                    gξ[i,j,ieq] = (a11*qx + a12*qy)*ωJ
                    gη[i,j,ieq] = (a21*qx + a22*qy)*ωJ
                end
            end
        end

        for j = 1:ngl, i = 1:ngl
            ip  = m.connijk[iel,i,j]
            ωJ  = m.ω[i]*m.ω[j]*m.Je[iel,i,j]
            a11 = m.dξdx[iel,i,j]; a12 = m.dξdy[iel,i,j]
            a21 = m.dηdx[iel,i,j]; a22 = m.dηdy[iel,i,j]
            for ieq = 1:neq
                dFdξ = dFdη = dGdξ = dGdη = 0.0
                for k = 1:ngl
                    dFdξ += m.dψ[k,i]*F[k,j,ieq]
                    dFdη += m.dψ[k,j]*F[i,k,ieq]
                    dGdξ += m.dψ[k,i]*F[k,j,neq+ieq]
                    dGdη += m.dψ[k,j]*F[i,k,neq+ieq]
                end
                r = -ωJ*((dFdξ*a11 + dFdη*a21) + (dGdξ*a12 + dGdη*a22) - S[i,j,ieq])
                if lvisc
                    tt = 0.0
                    for k = 1:ngl
                        tt += m.dψ[i,k]*gξ[k,j,ieq] + m.dψ[j,k]*gη[i,k,ieq]
                    end
                    r -= tt
                end
                RHS[ip,ieq] += r
            end
        end
    end

    for ieq = 1:neq, ip = 1:npoin
        RHS[ip,ieq] *= m.Minv[ip]
    end
    return RHS, uaux
end


#=================================================================================
                              WIRE THE CACHE BY HAND
=================================================================================#
function jacc_cache_from_test_mesh(m::TestMesh, qe::Matrix{Float64};
                                   neq = 3, lvisc = false, μ = Float64[],
                                   lsource = true, TW = Float64)
    T     = Float64                       # the host type, as in Jexpresso
    mixed = TW !== T
    p2e_ptr, p2e_e, p2e_i, p2e_j = _jacc_build_p2e(m.connijk, m.npoin, m.nelem, m.ngl)
    bn_ip, bn_ptr, bn_edge, bn_loc = _jacc_build_bnodes(m.poin_in_bdy_edge,
                                                        m.nedges_bdy, m.ngl)
    nb   = length(bn_ip)
    vdim = lvisc ? (m.nelem, m.ngl, m.ngl, neq) : (1,1,1,1)
    coef = lvisc ? collect(TW, μ) : Base.zeros(TW, neq)

    return St_jacc2d(npoin = m.npoin, nelem = m.nelem, ngl = m.ngl, neq = neq,
                     nbnode = nb,
                     TW = TW, mixed = mixed,
                     u_stage    = mixed ? Base.zeros(TW, m.npoin*neq) : Base.zeros(TW, 0),
                     bout       = _jacc_zeros(TW, max(nb,1), neq),
                     bout_host  = Base.zeros(TW, max(nb,1), neq),
                     bn_ip_host = bn_ip,
                     RHS_stage  = mixed ? Base.zeros(TW, m.npoin, neq) : Base.zeros(TW, 0, 0),
                     connijk = _jacc_copy(m.connijk),
                     coords = _jacc_copy(TW.(m.coords)),
                     x = _jacc_copy(TW.(m.x)), y = _jacc_copy(TW.(m.y)),
                     dξdx = _jacc_copy(TW.(m.dξdx)), dξdy = _jacc_copy(TW.(m.dξdy)),
                     dηdx = _jacc_copy(TW.(m.dηdx)), dηdy = _jacc_copy(TW.(m.dηdy)),
                     Je = _jacc_copy(TW.(m.Je)), dψ = _jacc_copy(TW.(m.dψ)),
                     ω = _jacc_copy(TW.(m.ω)), Minv = _jacc_copy(TW.(m.Minv)),
                     qe = _jacc_copy(TW.(qe)),
                     p2e_ptr = _jacc_copy(p2e_ptr), p2e_e = _jacc_copy(p2e_e),
                     p2e_i = _jacc_copy(p2e_i), p2e_j = _jacc_copy(p2e_j),
                     bn_ip = _jacc_copy(bn_ip), bn_ptr = _jacc_copy(bn_ptr),
                     bn_edge = _jacc_copy(bn_edge), bn_loc = _jacc_copy(bn_loc),
                     bdy_nx = _jacc_copy(TW.(m.nx)), bdy_ny = _jacc_copy(TW.(m.ny)),
                     u = _jacc_zeros(TW, m.npoin*neq),
                     uaux = _jacc_zeros(TW, m.npoin, neq+1),
                     flux = _jacc_zeros(TW, m.nelem, m.ngl, m.ngl, 2neq),
                     source = _jacc_zeros(TW, m.nelem, m.ngl, m.ngl, neq),
                     uprim = _jacc_zeros(TW, vdim...),
                     gξ = _jacc_zeros(TW, vdim...), gη = _jacc_zeros(TW, vdim...),
                     rhs_el = _jacc_zeros(TW, m.nelem, m.ngl, m.ngl, neq),
                     qbdy = _jacc_zeros(TW, max(nb,1), neq),
                     RHS = _jacc_zeros(TW, m.npoin, neq),
                     RHS_host = Base.zeros(T, m.npoin, neq),
                     visc_coeff = _jacc_copy(coef),
                     lvisc = lvisc, lsource = lsource, lpert = false,
                     PhysConst = nothing,
                     xmin = m.xmin, xmax = m.xmax, ymin = m.ymin, ymax = m.ymax)
end

# q (npoin × neq+1) → the flat ODE state the RHS is handed
function pack_u(q::Matrix{Float64}, npoin::Int, neq::Int)
    u = zeros(npoin*neq)
    for ieq = 1:neq, ip = 1:npoin
        u[(ieq-1)*npoin + ip] = q[ip,ieq]
    end
    return u
end


#=================================================================================
                                     MAIN
=================================================================================#
function main(nelx = 16, nely = 8, nop = 4)

    neq = 3
    μ   = [0.05, 0.05, 0.05]

    println()
    println("="^78)
    println(" JACC 2D shallow water — portable RHS test")
    println("="^78)
    @printf(" backend          : %s\n", JACC.backend)
    @printf(" julia threads    : %d\n", Threads.nthreads())
    @printf(" device array type: %s\n", nameof(typeof(JACC.zeros(Float64, 1))))

    m = build_test_mesh(nelx, nely, nop)
    @printf(" mesh             : %d x %d elements, nop = %d  →  npoin = %d, nelem = %d\n",
            nelx, nely, nop, m.npoin, m.nelem)
    println()

    q, qe = soliwave_state(m)
    nfail = 0

    #--- 1. well-balancedness: RHS(lake at rest) == 0 --------------------------
    let
        qrest = copy(qe)
        c  = jacc_cache_from_test_mesh(m, qe; neq = neq, lvisc = false)
        u  = pack_u(qrest, m.npoin, neq)
        jacc_rhs_2d!(c, u, 0.0)
        r  = maximum(abs, Array(c.RHS))
        ok = r < 1e-10
        nfail += ok ? 0 : 1
        @printf(" 1. lake at rest      max|RHS| = %.3e                    %s\n",
                r, ok ? "OK" : "FAIL")
    end

    #--- 2. JACC vs the serial reference, inviscid then viscous ---------------
    for (lv, tag) in ((false, "inviscid"), (true, "viscous "))
        c    = jacc_cache_from_test_mesh(m, qe; neq = neq, lvisc = lv, μ = μ)
        u    = pack_u(q, m.npoin, neq)
        jacc_rhs_2d!(c, u, 0.0)
        got  = Array(c.RHS)

        uaux = zeros(m.npoin, neq+1)
        for ieq = 1:neq, ip = 1:m.npoin
            uaux[ip,ieq] = q[ip,ieq]
        end
        want, _ = reference_rhs(m, uaux, qe, 0.0; neq = neq, lvisc = lv, μ = μ)

        scale = maximum(abs, want)
        aerr  = maximum(abs, got .- want)
        rerr  = aerr/max(scale, eps())
        # scale is the point of the line: a test that compared two zeros would
        # also report a zero difference, so print what the RHS actually is.
        ok    = rerr < 1e-12 && scale > 1e-6
        nfail += ok ? 0 : 1
        @printf(" 2. %s vs serial  max|RHS| = %.3e  max|Δ| = %.3e  rel = %.3e  %s\n",
                tag, scale, aerr, rerr, ok ? "OK" : "FAIL")
    end

    #--- 3. bitwise determinism ----------------------------------------------
    let
        c1 = jacc_cache_from_test_mesh(m, qe; neq = neq, lvisc = true, μ = μ)
        c2 = jacc_cache_from_test_mesh(m, qe; neq = neq, lvisc = true, μ = μ)
        jacc_rhs_2d!(c1, pack_u(q, m.npoin, neq), 0.0)
        jacc_rhs_2d!(c2, pack_u(q, m.npoin, neq), 0.0)
        ok = Array(c1.RHS) == Array(c2.RHS)
        nfail += ok ? 0 : 1
        @printf(" 3. determinism       two runs bitwise identical    %s\n", ok ? "OK" : "FAIL")
    end

    #--- 4. throughput --------------------------------------------------------
    let
        c = jacc_cache_from_test_mesh(m, qe; neq = neq, lvisc = true, μ = μ)
        u = pack_u(q, m.npoin, neq)
        jacc_rhs_2d!(c, u, 0.0)                      # warm up / compile
        nrep = 200
        t0 = time_ns()
        for _ = 1:nrep
            jacc_rhs_2d!(c, u, 0.0)
        end
        dt = (time_ns() - t0)/1e9
        @printf(" 4. throughput        %.1f RHS/s   (%.3f ms per call, %d reps)\n",
                nrep/dt, 1e3*dt/nrep, nrep)
        @printf("                      %.2f Mdof/s\n", nrep*m.npoin*neq/dt/1e6)
    end

    println()
    if nfail == 0
        println(" ALL CHECKS PASSED")
    else
        println(" ", nfail, " CHECK(S) FAILED")
    end
    println("="^78)
    println()
    return nfail
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    args = length(ARGS) >= 3 ? parse.(Int, ARGS[1:3]) : [16, 8, 4]
    exit(main(args...) == 0 ? 0 : 1)
end
