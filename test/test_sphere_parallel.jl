#---------------------------------------------------------------------------------
# test/test_sphere_parallel.jl — the spherical shell under MPI.
#
#   mpiexec -n 4 julia --project=. test/test_sphere_parallel.jl
#
# (any rank count works, including 1; from a Julia session,
#  `using MPI; run(`$(mpiexec()) -n 4 $(Base.julia_cmd()[1]) --project=. test/test_sphere_parallel.jl`)`.)
#
# WHAT BREAKS WHEN A CONTINUOUS-GALERKIN MANIFOLD SOLVER IS PARTITIONED, and
# what each test below pins down.
#
# Every rank integrates only the elements it owns. A node on a partition seam is
# touched by elements on BOTH sides, so each rank holds a PARTIAL sum of every
# assembled quantity there — the mass matrix, the RHS, the filter average, the
# vorticity. Completing those sums across ranks (assemble_mpi!, the cross-rank
# half of the direct stiffness summation) is the whole of what makes the shell
# parallel, and getting it wrong does NOT crash: it leaves a seam of nodes
# carrying a fraction of their divergence, i.e. a ring of spurious forcing
# threaded through the mesh, which looks like a plausible solution.
#
# So the tests are the invariants that a broken seam violates:
#
#   1. THE AREA. Σ M[ip] over OWNED nodes must be 4πR². It fails high if the
#      ownership test is missing (shared nodes counted once per rank) and low if
#      the assembly is missing (shared nodes carry only their local share).
#
#   2. AGREEMENT ON SHARED NODES. A node held by k ranks must carry the SAME
#      mass and the SAME RHS on all k of them; after assemble_mpi! the owner's
#      total is sent back to the mirrors, so this is exact, not approximate.
#      Checked by gathering (global id, value) pairs on rank 0.
#
#   3. ∫ ∂φ/∂t dΩ = 0. The shell is closed, so the weak divergence of the
#      continuity flux integrates to zero over it: Σ M[ip]·RHS[ip,1] must vanish
#      to round-off. THIS is the test that a partial sum at a seam cannot pass —
#      the missing element contributions do not cancel — and it needs no serial
#      reference run to compare against.
#
#   4. The seven metric checks of check_sphere_metrics, which must return the
#      same verdict on every rank (they reduce before they judge).
#
# The fixture is the .msh that ships with the SWsphere case, so the test needs
# no gmsh binary.
#---------------------------------------------------------------------------------

using Test
using Jexpresso
using Jexpresso: mod_mesh_mesh_driver, build_sphere_metrics, check_sphere_metrics,
                 build_sphere_params, sphere_rhs!, sphere_filter!,
                 sphere_relative_vorticity!
using PartitionedArrays, MPI
using Printf

const CASE_MSH = joinpath(@__DIR__, "..", "problems", "ShallowWater", "SWsphere",
                          "cubed_sphere.msh")

# mod_inputs_user_inputs! reads three module globals that run.jl normally sets
# from the command line before it includes a case.
@eval Jexpresso begin
    parsed_equations           = "ShallowWater"
    parsed_equations_case_name = "SWsphere"
    user_input_file            = "test/test_sphere_parallel.jl"
end

# The case's hooks — user_flux!, user_source!, initialize, user_uout! — are
# ordinary functions that run.jl `include`s INTO the Jexpresso module. This test
# calls the RHS, so it needs them; do what run.jl does.
const CASE_DIR = joinpath(@__DIR__, "..", "problems", "ShallowWater", "SWsphere")
for _f in ("user_flux.jl", "user_source.jl", "user_bc.jl",
           "initialize.jl", "user_primitives.jl")
    Base.include(Jexpresso, joinpath(CASE_DIR, _f))
end

function shell_inputs(nop)
    inputs = Dict{Symbol,Any}(
        :lread_gmsh           => true,
        :gmsh_filename        => CASE_MSH,
        :nop                  => nop,
        :interpolation_nodes  => "lgl",
        :backend              => Jexpresso.CPU(),
        :lspherical_shell     => true,
        :sphere_radius        => 6.37122e6,
        :lproject_to_sphere   => true,
        :lgrid_only           => true,
        :sphere_metrics       => :radial,
    )
    Jexpresso.mod_inputs_user_inputs!(inputs, 1)   # rank 1 => no banner spam
    inputs[:gmsh_filename]  = CASE_MSH
    inputs[:nop]            = nop
    inputs[:sphere_metrics] = :radial
    return inputs
end

#
# How many DISTINCT global nodes the ranks hold between them. Not the same as
# mesh.gnpoin, which is max(ip2gip) — the high-order global numbering leaves
# gaps, so on the shipped grid at nop = 4 the ids run to 9762 while there are
# 9602 nodes. The count below is the one that has to match "every node is owned
# by exactly one rank".
#
function count_distinct_gips(gips, comm)
    counts = MPI.Allgather(Int32(length(gips)), comm)
    allg   = MPI.Gatherv!(collect(Int64, gips), MPI.VBuffer(
                          Vector{Int64}(undef, sum(counts)), counts), 0, comm)
    n = MPI.Comm_rank(comm) == 0 ? length(Set(allg)) : 0
    return MPI.bcast(n, 0, comm)
end

#
# Gather (global node id, value) from every rank and report the largest
# disagreement between two ranks that hold the same node — 0.0 when nothing is
# shared, which is what a one-rank run sees.
#
function worst_shared_disagreement(gips, vals, comm)
    counts = MPI.Allgather(Int32(length(gips)), comm)
    allg   = MPI.Gatherv!(collect(Int64, gips), MPI.VBuffer(
                          Vector{Int64}(undef, sum(counts)), counts), 0, comm)
    allv   = MPI.Gatherv!(collect(Float64, vals), MPI.VBuffer(
                          Vector{Float64}(undef, sum(counts)), counts), 0, comm)
    worst  = 0.0
    if MPI.Comm_rank(comm) == 0
        seen = Dict{Int64,Float64}()
        for k in eachindex(allg)
            g = allg[k]; v = allv[k]
            if haskey(seen, g)
                s = max(abs(v), abs(seen[g]))
                worst = max(worst, s > 0 ? abs(v - seen[g])/s : abs(v - seen[g]))
            else
                seen[g] = v
            end
        end
    end
    return MPI.bcast(worst, 0, comm)
end


with_mpi() do distribute

    comm   = MPI.COMM_WORLD
    rank   = MPI.Comm_rank(comm)
    nparts = MPI.Comm_size(comm)
    isfile(CASE_MSH) || error("missing $CASE_MSH — it ships with the SWsphere case")

    rank == 0 && @info "spherical shell on $nparts MPI rank(s)"

    @testset "spherical shell on $nparts rank(s)" begin

        nop     = 4
        inputs  = shell_inputs(nop)
        mesh, _ = mod_mesh_mesh_driver(inputs, nparts, distribute)
        @test mesh.lmanifold

        metrics = build_sphere_metrics(mesh, inputs; verbose = false)

        #--- 4. the seven metric checks, same verdict everywhere
        ok = check_sphere_metrics(mesh, metrics; verbose = (rank == 0))
        @test ok
        @test MPI.Allreduce(ok ? 1 : 0, MPI.MIN, comm) == 1   # nobody disagreed

        npoin = Int(mesh.npoin)
        owned = [ip for ip = 1:npoin if mesh.gip2owner[ip] == rank]

        #--- 1. the area, over owned nodes only
        area = MPI.Allreduce(sum(@view metrics.M[owned]), MPI.SUM, comm)
        aref = 4π*mesh.radius^2
        rank == 0 && @printf(" #   area   Σ M[owned] = %.10e , 4πR² = %.10e , rel err %.2e\n",
                             area, aref, abs(area - aref)/aref)
        @test abs(area - aref)/aref < 1.0e-10

        # OWNERSHIP IS A PARTITION of the node set: the owned counts add up to
        # the number of distinct nodes, so no node is owned twice (which would
        # double-count it in every integral) and none is orphaned (which would
        # drop it from every integral).
        nown_total = MPI.Allreduce(length(owned), MPI.SUM, comm)
        ndistinct  = count_distinct_gips(@view(mesh.ip2gip[1:npoin]), comm)
        # Reduced FIRST, printed second. An MPI call inside a `rank == 0 &&`
        # is evaluated on rank 0 only, and hangs every other rank.
        nlocal_total = MPI.Allreduce(npoin, MPI.SUM, comm)
        rank == 0 && @printf(" #   nodes  owned Σ = %d , distinct global = %d , local Σ = %d\n",
                             nown_total, ndistinct, nlocal_total)
        @test nown_total == ndistinct

        #--- 2. shared nodes agree on the assembled mass
        dM = worst_shared_disagreement(@view(mesh.ip2gip[1:npoin]),
                                       @view(metrics.M[1:npoin]), comm)
        rank == 0 && @printf(" #   mass matrix: worst disagreement on a shared node = %.2e\n", dM)
        @test dM < 1.0e-13

        #--- the RHS of the shipped initial condition
        q  = Jexpresso.initialize(mesh.SD, 0, mesh, inputs, mktempdir(), Float64)
        sp = build_sphere_params(mesh, metrics, inputs; neqs = q.neqs)
        RHS = zeros(Float64, npoin, q.neqs + 1)
        sphere_rhs!(RHS, q.qn, q.qe, mesh, metrics, sp, inputs[:SOL_VARS_TYPE])

        #--- 2. (again) shared nodes agree on the assembled RHS
        dR = worst_shared_disagreement(@view(mesh.ip2gip[1:npoin]),
                                       @view(RHS[1:npoin, 1]), comm)
        rank == 0 && @printf(" #   RHS[φ]:      worst disagreement on a shared node = %.2e\n", dR)
        @test dR < 1.0e-12

        #--- 3. ∫ ∂φ/∂t dΩ = 0 on a closed shell. The seam test.
        dmass = MPI.Allreduce(sum(metrics.M[ip]*RHS[ip,1] for ip in owned), MPI.SUM, comm)
        scale = MPI.Allreduce(sum(metrics.M[ip]*abs(RHS[ip,1]) for ip in owned), MPI.SUM, comm)
        rank == 0 && @printf(" #   ∫∂φ/∂t dΩ = %.3e against Σ M|∂φ/∂t| = %.3e  (ratio %.2e)\n",
                             dmass, scale, abs(dmass)/scale)
        @test abs(dmass)/scale < 1.0e-12

        #--- the filter and the vorticity go through the same assembly, so give
        #    them the same shared-node treatment
        ζ = zeros(Float64, npoin)
        sphere_relative_vorticity!(ζ, q.qn, mesh, metrics, sp)
        dζ = worst_shared_disagreement(@view(mesh.ip2gip[1:npoin]), ζ, comm)
        rank == 0 && @printf(" #   vorticity:   worst disagreement on a shared node = %.2e\n", dζ)
        @test dζ < 1.0e-12

        if sp.lfilter
            qf = copy(q.qn)
            sphere_filter!(qf, mesh, metrics, sp)
            dF = worst_shared_disagreement(@view(mesh.ip2gip[1:npoin]), @view(qf[1:npoin,1]), comm)
            rank == 0 && @printf(" #   filtered φ:  worst disagreement on a shared node = %.2e\n", dF)
            @test dF < 1.0e-13
        end
    end
end
