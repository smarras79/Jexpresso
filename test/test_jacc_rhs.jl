#---------------------------------------------------------------------------------
# test/test_jacc_rhs.jl — the portable JACC.jl right-hand side
#                         (src/kernel/operators/rhs_jacc.jl).
#
#   julia --project=. -t 4 test/test_jacc_rhs.jl
#
# Deliberately does NOT `using Jexpresso`. The kernels in rhs_jacc.jl take plain
# arrays and the four device callbacks of the case; nothing in them needs the
# module, its 60-package dependency tree, MPI, or a gmsh grid. Keeping the test
# free of all that is what makes it runnable in seconds on a machine where a JACC
# backend is being brought up for the first time — which is exactly when it is
# wanted. The mesh, the reference implementation and the SoliWaveIsland callbacks
# come from problems/JACCtestz/JACC_2d_swe.jl, which this file drives.
#
# WHAT IS ASSERTED
#
#   * the discrete operator annihilates the lake at rest — the well-balanced
#     flux/source split of problems/ShallowWater/SoliWaveIsland/user_flux.jl is a
#     property of the DISCRETISATION, not just of the continuous equations, and
#     this is where it gets checked;
#   * the JACC RHS equals an independently written serial evaluation of the same
#     weak form, inviscid and viscous, to the last bit;
#   * two evaluations agree bitwise. They must: the direct stiffness summation is
#     a gather over a fixed CSR map, not an atomic scatter, so the summation
#     order does not depend on how the threads (or the warps) happened to be
#     scheduled. This is the property that makes a GPU run comparable with a CPU
#     run at all, and it is the reason the kernels do not use atomics;
#   * the point→element map and the boundary-node map really are inverses of
#     connijk / poin_in_bdy_edge;
#   * the free-slip wall leaves no normal momentum behind, at the corners too —
#     the case a per-edge parallelisation gets wrong.
#
# The mesh is curvilinear (warped interior, straight walls), so dξdy and dηdx are
# non-zero and a transposed metric or a swapped index cannot pass.
#---------------------------------------------------------------------------------

using Test
using JACC

include(joinpath(@__DIR__, "..", "problems", "JACCtestz", "JACC_2d_swe.jl"))

const NEQ = 3
const MU  = [0.05, 0.05, 0.05]

@testset verbose = true "JACC portable RHS (2D shallow water)" begin

    m      = build_test_mesh(12, 6, 4)
    q, qe  = soliwave_state(m)

    @testset "backend in use" begin
        # Not an assertion about WHICH backend — it is a log line, so a CI run on
        # a GPU box says so in its output instead of silently testing threads.
        @info "JACC backend" backend=JACC.backend threads=Threads.nthreads() npoin=m.npoin nelem=m.nelem
        @test m.npoin > 0
    end

    @testset "connectivity maps are inverses" begin
        ptr, pe, pi_, pj = _jacc_build_p2e(m.connijk, m.npoin, m.nelem, m.ngl)
        @test length(ptr) == m.npoin + 1
        @test ptr[end] - 1 == m.nelem*m.ngl*m.ngl
        # every listed (iel,i,j) must map back to the node that owns the slot
        ok = true
        for ip = 1:m.npoin, mm = ptr[ip]:(ptr[ip+1]-1)
            m.connijk[pe[mm], pi_[mm], pj[mm]] == ip || (ok = false)
        end
        @test ok

        bn_ip, bn_ptr, bn_edge, bn_loc = _jacc_build_bnodes(m.poin_in_bdy_edge,
                                                            m.nedges_bdy, m.ngl)
        @test length(bn_ip) == length(unique(vec(m.poin_in_bdy_edge)))
        @test bn_ptr[end] - 1 == length(m.poin_in_bdy_edge)   # every slot accounted for
        okb = true
        for b = 1:length(bn_ip), mm = bn_ptr[b]:(bn_ptr[b+1]-1)
            m.poin_in_bdy_edge[bn_edge[mm], bn_loc[mm]] == bn_ip[b] || (okb = false)
        end
        @test okb
        # a corner of the rectangle belongs to two edges — the whole reason the
        # boundary kernel is parallelised over NODES and not over (edge, node)
        @test maximum(diff(bn_ptr)) == 2
    end

    @testset "lake at rest is an exact steady state" begin
        c = jacc_cache_from_test_mesh(m, qe; neq = NEQ, lvisc = false)
        jacc_rhs_2d!(c, pack_u(copy(qe), m.npoin, NEQ), 0.0)
        @test maximum(abs, Array(c.RHS)) < 1e-12
    end

    @testset "matches the serial reference" begin
        uaux = zeros(m.npoin, NEQ+1)
        for ieq = 1:NEQ, ip = 1:m.npoin
            uaux[ip,ieq] = q[ip,ieq]
        end

        for lvisc in (false, true)
            c = jacc_cache_from_test_mesh(m, qe; neq = NEQ, lvisc = lvisc, μ = MU)
            jacc_rhs_2d!(c, pack_u(q, m.npoin, NEQ), 0.0)
            got  = Array(c.RHS)
            want, _ = reference_rhs(m, uaux, qe, 0.0; neq = NEQ, lvisc = lvisc, μ = MU)

            # the state is a solitary wave running onto a cone: the RHS had better
            # not be zero, or "they agree" would mean nothing
            @test maximum(abs, want) > 1e-3
            @test got == want
        end
    end

    @testset "mixed precision (Float32 residual)" begin
        # This is the mode Metal needs, and it is exercised here on the threads
        # backend so it is covered on machines that have no Apple GPU. Nothing in
        # the staging is backend-specific: the cache carries TW-typed device
        # arrays and the driver converts on the two host hops, whatever TW is.

        c32 = jacc_cache_from_test_mesh(m, qe; neq = NEQ, lvisc = true, μ = MU,
                                        TW = Float32)
        @test c32.mixed
        @test eltype(Array(c32.RHS)) === Float32
        @test eltype(c32.RHS_host)   === Float64   # the host mirror stays double

        # 1. the lake at rest is STILL an exact steady state. H and He are equal
        #    before the cast so they are equal after it, and the well-balanced
        #    flux/source split cancels term by term in any precision. If this ever
        #    fails, the conversion is being applied to H and He inconsistently.
        crest = jacc_cache_from_test_mesh(m, qe; neq = NEQ, lvisc = false,
                                          TW = Float32)
        jacc_rhs_2d!(crest, pack_u(copy(qe), m.npoin, NEQ), 0.0)
        @test maximum(abs, Array(crest.RHS)) == 0.0

        # 2. against the Float64 residual, single-precision accuracy — and no
        #    worse. The bound is loose on purpose: what is being asserted is that
        #    the answer degrades by eps32 and not by something structural.
        c64 = jacc_cache_from_test_mesh(m, qe; neq = NEQ, lvisc = true, μ = MU)
        jacc_rhs_2d!(c64, pack_u(q, m.npoin, NEQ), 0.0)
        jacc_rhs_2d!(c32, pack_u(q, m.npoin, NEQ), 0.0)
        r64 = Array(c64.RHS)
        r32 = Float64.(Array(c32.RHS))
        rel = maximum(abs, r32 .- r64) / maximum(abs, r64)
        @test rel < 1e-5
        @test rel > 0            # it really did run in single, not silently in double

        # 3. the state stays Float64 and only the boundary nodes are written back
        u  = pack_u(q, m.npoin, NEQ)
        b4 = copy(u)
        jacc_rhs_2d!(c32, u, 0.0)
        @test eltype(u) === Float64
        bnodes  = Set(vec(m.poin_in_bdy_edge))
        changed = findall(!=(0.0), u .- b4)
        @test !isempty(changed)
        @test all(k -> ((k - 1) % m.npoin + 1) in bnodes, changed)

        # 4. still deterministic
        c32b = jacc_cache_from_test_mesh(m, qe; neq = NEQ, lvisc = true, μ = MU,
                                         TW = Float32)
        jacc_rhs_2d!(c32b, pack_u(q, m.npoin, NEQ), 0.0)
        @test Array(c32.RHS) == Array(c32b.RHS)
    end

    @testset "bitwise reproducible" begin
        r = map(1:2) do _
            c = jacc_cache_from_test_mesh(m, qe; neq = NEQ, lvisc = true, μ = MU)
            jacc_rhs_2d!(c, pack_u(q, m.npoin, NEQ), 0.0)
            Array(c.RHS)
        end
        @test r[1] == r[2]
    end

    @testset "free-slip wall, corners included" begin
        # The boundary values are applied to `uaux` and then mirrored back into the
        # state vector, which is what uaux2u! does at the tail of the CPU's
        # Dirichlet builder. On these straight walls n is (±1,0) or (0,±1), so
        # no-penetration means the wall-normal momentum is zero — at a corner both
        # components are.
        c = jacc_cache_from_test_mesh(m, qe; neq = NEQ, lvisc = false)
        u = pack_u(q, m.npoin, NEQ)
        before = copy(u)
        jacc_rhs_2d!(c, u, 0.0)

        # the corrected field really did reach `u`, and only at the wall
        @test u != before
        changed = findall(!=(0.0), u .- before)
        bnodes  = Set(vec(m.poin_in_bdy_edge))
        @test all(ip -> ((ip - 1) % m.npoin + 1) in bnodes, changed)

        ua = Array(c.uaux)
        worst = 0.0
        for iedge = 1:m.nedges_bdy, k = 1:m.ngl
            ip = m.poin_in_bdy_edge[iedge,k]
            worst = max(worst, abs(m.nx[iedge,k]*ua[ip,2] + m.ny[iedge,k]*ua[ip,3]))
        end
        # 1e-6 and not 1e-14: the CPU path declines any boundary correction the
        # AlmostEqual gate calls negligible, and this path reproduces that exactly.
        @test worst < 2e-6

        npx     = m.nelx*(m.ngl - 1) + 1
        corners = [1, npx]                          # the two south-wall corners
        for ip in corners
            @test abs(ua[ip,2]) < 2e-6
            @test abs(ua[ip,3]) < 2e-6
            # and the same values are what the integrator now holds
            @test u[1*m.npoin + ip] == ua[ip,2]
            @test u[2*m.npoin + ip] == ua[ip,3]
        end
    end
end
