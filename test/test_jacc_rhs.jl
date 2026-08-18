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

    @testset "bitwise reproducible" begin
        r = map(1:2) do _
            c = jacc_cache_from_test_mesh(m, qe; neq = NEQ, lvisc = true, μ = MU)
            jacc_rhs_2d!(c, pack_u(q, m.npoin, NEQ), 0.0)
            Array(c.RHS)
        end
        @test r[1] == r[2]
    end

    @testset "free-slip wall, corners included" begin
        # The boundary kernel writes back into the state vector it is handed, so
        # `u` after the call carries the Dirichlet values. On these straight walls
        # n is (±1,0) or (0,±1), so no-penetration means the wall-normal momentum
        # component is exactly zero — and at a corner BOTH components are.
        c = jacc_cache_from_test_mesh(m, qe; neq = NEQ, lvisc = false)
        u = pack_u(q, m.npoin, NEQ)
        jacc_rhs_2d!(c, u, 0.0)

        worst = 0.0
        for iedge = 1:m.nedges_bdy, k = 1:m.ngl
            ip = m.poin_in_bdy_edge[iedge,k]
            Hu = u[1*m.npoin + ip]      # q = [H, Hu, Hv] ⇒ momentum starts at neqs 2
            Hv = u[2*m.npoin + ip]
            worst = max(worst, abs(m.nx[iedge,k]*Hu + m.ny[iedge,k]*Hv))
        end
        @test worst < 1e-14

        # a corner sees two normals, so both components must have been zeroed
        npx     = m.nelx*(m.ngl - 1) + 1
        corners = [1, npx]                          # the two south-wall corners
        for ip in corners
            @test abs(u[1*m.npoin + ip]) < 1e-14
            @test abs(u[2*m.npoin + ip]) < 1e-14
        end
    end
end
