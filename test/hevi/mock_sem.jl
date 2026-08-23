#=============================================================================
 mock_sem.jl -- a minimal spectral-element stand-in for the HEVI unit tests.

 The HEVI sources touch a small, well-defined slice of Jexpresso: a mesh with
 `connijk`/`coords`/`ip2gip`, the metric terms, the LGL basis and weights, the
 lumped inverse mass matrix, the DSS routines and the MPI assembler. This file
 provides exactly that slice for a uniform Cartesian box, so the column
 topology, the MPI redistribution, the probed band assembly and the stage
 solve can be tested for real -- on several ranks, with columns deliberately
 split across them -- without building the whole dependency stack.

 It is a test fixture, not a second implementation: `DSS_rhs!` and
 `divide_by_mass_matrix!` are copies of the production versions, and the
 assembler has the same contract (sum contributions at shared nodes, then
 make every copy equal).
=============================================================================#

using MPI, LinearAlgebra, Printf

struct NSD_1D end
struct NSD_3D end
struct ContGal end

Base.@kwdef struct PhysicalConst{T}
    Rair::T = 287.0
    cp::T   = 1004.0
    cv::T   = 718.0
    γ::T    = T(1004.0/718.0)
    pref::T = 101200.0
    C0::T   = T((287.0^(1004.0/718.0)) / 101200.0^((1004.0/718.0) - 1.0))
    g::T    = 9.80616
    μ_mol::T = 1.8e-5
    Pr_t::T  = 0.7
end
@inline perfectGasLaw_ρθtoP(P::PhysicalConst, ρ::Real, θ::Real) = P.C0 * (ρ*θ)^P.γ

#-----------------------------------------------------------------------------
# Legendre-Gauss-Lobatto nodes, weights and differentiation matrix
#-----------------------------------------------------------------------------
function legendre(n::Int, x::Float64)
    n == 0 && return (1.0, 0.0)
    n == 1 && return (x, 1.0)
    pm2, pm1 = 1.0, x
    dm2, dm1 = 0.0, 1.0
    p = pm1; d = dm1
    for k = 2:n
        p = ((2k-1)*x*pm1 - (k-1)*pm2) / k
        d = dm2 + (2k-1)*pm1
        pm2, pm1 = pm1, p
        dm2, dm1 = dm1, d
    end
    return (p, d)
end

"""
    lgl(ngl) -> (ξ, ω, dψ)

`dψ[m,k]` is dℓ_m/dξ evaluated at node k, matching Jexpresso's convention.
"""
function lgl(ngl::Int)
    N = ngl - 1
    ξ = zeros(ngl); ξ[1] = -1.0; ξ[end] = 1.0
    for i = 2:ngl-1                      # Chebyshev start, Newton on P'_N
        x = -cos(pi * (i-1) / N)
        for _ = 1:100
            p, d = legendre(N, x)
            # P'_N = d ; need root of P'_N, so Newton on d with d' from the ODE
            dd = (2x*d - N*(N+1)*p) / (1 - x^2)
            dx = -d / dd
            x += dx
            abs(dx) < 1e-15 && break
        end
        ξ[i] = x
    end
    ω = [2.0 / (N*(N+1) * legendre(N, ξ[i])[1]^2) for i = 1:ngl]
    D = zeros(ngl, ngl)                  # D[k,m] = dℓ_m/dξ at ξ_k
    for k = 1:ngl, m = 1:ngl
        if k != m
            D[k,m] = legendre(N, ξ[k])[1] / (legendre(N, ξ[m])[1] * (ξ[k] - ξ[m]))
        end
    end
    D[1,1]       = -N*(N+1)/4
    D[ngl,ngl]   =  N*(N+1)/4
    dψ = permutedims(D)                  # dψ[m,k]
    return ξ, ω, dψ
end

#-----------------------------------------------------------------------------
# The MPI assembler: same contract as Jexpresso's, implemented the slow,
# obviously-correct way (a global reduction) because this is a test.
#-----------------------------------------------------------------------------
struct MockAssembler
    comm::MPI.Comm
    gnpoin::Int
    ip2gip::Vector{Int}
    G::Matrix{Float64}          # preallocated: this sits on the RHS hot path
end
function setup_assembler(SD, a, ip2gip, gip2owner)
    gn = MPI.Allreduce(maximum(ip2gip), MPI.MAX, MPI.COMM_WORLD)
    MockAssembler(MPI.COMM_WORLD, gn, copy(ip2gip), zeros(Float64, gn, size(a, 2)))
end

function assemble_mpi!(a::AbstractMatrix, c::MockAssembler)
    m = size(a, 2)
    G = m == size(c.G, 2) ? fill!(c.G, 0.0) : zeros(Float64, c.gnpoin, m)
    @inbounds for (ip, g) in enumerate(c.ip2gip), j = 1:m
        G[g, j] += a[ip, j]
    end
    MPI.Allreduce!(G, MPI.SUM, c.comm)
    @inbounds for (ip, g) in enumerate(c.ip2gip), j = 1:m
        a[ip, j] = G[g, j]
    end
    return a
end

# --- production copies -------------------------------------------------------
function DSS_rhs!(RHS, rhs_el, connijk, nelem, ngl, neqs, ::NSD_3D, ::ContGal)
    for ieq = 1:neqs, iel = 1:nelem, k = 1:ngl, j = 1:ngl, i = 1:ngl
        RHS[connijk[iel,i,j,k], ieq] += rhs_el[iel,i,j,k,ieq]
    end
end
function divide_by_mass_matrix!(RHS, RHSaux, Minv::AbstractVector, neqs, npoin, ::ContGal)
    for ip = 1:npoin
        RHS[ip] = Minv[ip] * RHS[ip]
    end
end

#-----------------------------------------------------------------------------
# The synthetic mesh
#-----------------------------------------------------------------------------
struct MockParts; comm::MPI.Comm; end
struct MockMesh
    npoin::Int; nelem::Int; ngl::Int
    connijk::Array{Int,4}
    coords::Matrix{Float64}
    xmin::Float64; xmax::Float64; ymin::Float64; ymax::Float64; zmin::Float64; zmax::Float64
    ip2gip::Vector{Int}; gip2owner::Vector{Int}
    SD::NSD_3D; parts::MockParts
end
struct MockMetrics
    Je::Array{Float64,4}
    dξdx::Array{Float64,4}; dηdy::Array{Float64,4}
    dζdx::Array{Float64,4}; dζdy::Array{Float64,4}; dζdz::Array{Float64,4}
end

"""
    UntypedMetrics

The same metrics with NO field type annotations -- i.e. every field `::Any`.

This is not a straw man: `St_metrics` and `St_mesh` in Jexpresso are
`Base.@kwdef mutable struct`s written exactly this way, so `metrics.dζdz` and
`mesh.connijk` really are `::Any` in a production run. Reading them inside an
element loop makes every index a boxed dynamic dispatch, and the operator ran
60x slower against the real structs than against the typed mocks here --
invisible to every test in this directory until it showed up as HEVI being
slower than the explicit scheme on a real case.

`build_mock_params(; untyped_metrics = true)` reproduces that, so the typed
function barrier in hevi_apply_A! has something to be tested against.
"""
mutable struct UntypedMetrics
    Je
    dξdx; dηdy
    dζdx; dζdy; dζdz
end
struct MockBasis; dψ::Matrix{Float64}; end
struct MockQP; qe::Matrix{Float64}; end

"""
    build_mock_params(; nelx, nely, nelz, p, Lx, Ly, Lz, comm) -> params

Uniform Cartesian mesh of `nelx x nely x nelz` order-`p` elements, partitioned
by splitting the element list into contiguous per-rank chunks.

Elements are numbered with the vertical index fastest, so a contiguous chunk
is a run of whole columns plus a partial one at each end -- deliberately the
same shape a p4est space-filling curve produces over a column-ordered coarse
mesh, and the case the redistribution has to get right.
"""
function build_mock_params(; nelx = 2, nely = 2, nelz = 5, p = 2,
                             Lx = 400.0, Ly = 400.0, Lz = 1000.0,
                             comm = MPI.COMM_WORLD, θ0 = 300.0, ρ0 = 1.2,
                             base = :exp, untyped_metrics = false)

    rank   = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)
    ngl    = p + 1
    ξ, ω, dψ = lgl(ngl)

    npx, npy, npz = nelx*p + 1, nely*p + 1, nelz*p + 1
    gnpoin = npx * npy * npz
    gid(ix, iy, iz) = (iz-1)*npx*npy + (iy-1)*npx + ix

    hx, hy, hz = Lx/nelx, Ly/nely, Lz/nelz

    # element list, vertical fastest
    els = [(ex, ey, ez) for ex = 1:nelx for ey = 1:nely for ez = 1:nelz]
    nelG = length(els)
    lo = div(nelG * rank, nranks) + 1
    hi = div(nelG * (rank+1), nranks)
    mine = els[lo:hi]
    nelem = length(mine)

    # local numbering of the gids this rank touches
    g2l = Dict{Int,Int}(); ip2gip = Int[]
    conn = zeros(Int, nelem, ngl, ngl, ngl)
    for (iel, (ex, ey, ez)) in enumerate(mine)
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            g = gid((ex-1)*p + i, (ey-1)*p + j, (ez-1)*p + k)
            l = get(g2l, g, 0)
            if l == 0
                push!(ip2gip, g); l = length(ip2gip); g2l[g] = l
            end
            conn[iel, i, j, k] = l
        end
    end
    npoin = length(ip2gip)

    coords = zeros(3, npoin)
    for (iel, (ex, ey, ez)) in enumerate(mine)
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = conn[iel, i, j, k]
            coords[1, ip] = (ex-1)*hx + 0.5*hx*(1 + ξ[i])
            coords[2, ip] = (ey-1)*hy + 0.5*hy*(1 + ξ[j])
            coords[3, ip] = (ez-1)*hz + 0.5*hz*(1 + ξ[k])
        end
    end

    # lowest rank holding a gid owns it
    ownG = fill(typemax(Int32), gnpoin)
    for g in ip2gip; ownG[g] = Int32(rank); end
    MPI.Allreduce!(ownG, MPI.MIN, comm)
    gip2owner = [Int(ownG[g]) for g in ip2gip]

    Je   = fill(0.125*hx*hy*hz, nelem, ngl, ngl, ngl)
    zer  = zeros(nelem, ngl, ngl, ngl)
    _mtype  = untyped_metrics ? UntypedMetrics : MockMetrics
    metrics = _mtype(Je,
                     fill(2.0/hx, nelem, ngl, ngl, ngl),
                     fill(2.0/hy, nelem, ngl, ngl, ngl),
                     copy(zer), copy(zer), fill(2.0/hz, nelem, ngl, ngl, ngl))

    # lumped mass, assembled across ranks
    M = zeros(npoin, 1)
    for iel = 1:nelem, k = 1:ngl, j = 1:ngl, i = 1:ngl
        M[conn[iel,i,j,k], 1] += ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
    end
    mesh = MockMesh(npoin, nelem, ngl, conn, coords, 0.0, Lx, 0.0, Ly, 0.0, Lz,
                    ip2gip, gip2owner, NSD_3D(), MockParts(comm))
    asm = setup_assembler(NSD_3D(), M, ip2gip, gip2owner)
    assemble_mpi!(M, asm)
    Minv = [1.0 / M[ip, 1] for ip = 1:npoin]

    # Reference state.
    #
    #   :exp     a cheap exponential density profile at constant θ. NOT
    #            hydrostatic, which is fine for the linear-operator tests --
    #            they only need a smooth state to linearise about.
    #   :neutral the exact neutrally stratified hydrostatic profile, which is
    #            what a rising-thermal-bubble benchmark is posed on. Getting
    #            this right matters there: the whole point of a warm bubble is
    #            that the perturbation is the only thing driving motion, so any
    #            imbalance in the base state shows up as spurious rising.
    PC = PhysicalConst{Float64}()
    qe = zeros(npoin, 5)
    for ip = 1:npoin
        z = coords[3, ip]
        if base === :neutral
            pbar = PC.pref * (1.0 - PC.g*z/(PC.cp*θ0))^(PC.cp/PC.Rair)
            ρ    = (1.0/θ0) * (pbar/PC.C0)^(1.0/PC.γ)
        else
            ρ = ρ0 * exp(-z / 8500.0)
        end
        qe[ip, 1] = ρ
        qe[ip, 5] = ρ * θ0
    end

    return (mesh = mesh, metrics = metrics, basis = MockBasis(dψ), ω = ω,
            Minv = Minv, SD = NSD_3D(), AD = ContGal(), qp = MockQP(qe), neqs = 5,
            inputs = Dict{Symbol,Any}(), hx = hx, hy = hy, hz = hz, ngl = ngl)
end
