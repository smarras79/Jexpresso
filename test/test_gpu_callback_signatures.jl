#---------------------------------------------------------------------------------
# test/test_gpu_callback_signatures.jl — every case's device callbacks must have
# the arity the kernels actually call them with.
#
#   julia --project=. test/test_gpu_callback_signatures.jl
#
# NO `using Jexpresso`: this only parses files, so it runs in a bare Julia in
# under a second.
#
# WHAT IT GUARDS, AND WHY IT EXISTS
#
# `user_bc_dirichlet_gpu` and friends are resolved at the CALL SITE, inside a
# kernel, at run time. Nothing links them at load time, nothing type-checks them,
# and no case exercises them unless someone actually runs that case on a GPU. So
# a signature that does not match the kernel is not a compile error, not a test
# failure, and not visible in review — it is a MethodError that fires the first
# time somebody points a GPU at that deck, possibly years later.
#
# That is not hypothetical. When the coordinate arguments were migrated from
# mesh.x/mesh.y/mesh.z to mesh.coords, the kernels kept passing `x[ip], y[ip]`
# while 46 of the 70 case directories had already been rewritten to take a single
# `coords` slot. Every one of those 46 cases had a GPU boundary path that could
# only ever throw. It stayed that way because there was nothing like this file.
#
# The check is deliberately crude — count the parameters, compare against what
# the kernel passes. It cannot tell you the body is right. It can tell you the
# call cannot possibly work, which is the failure mode that actually happened.
#---------------------------------------------------------------------------------

using Test

const ROOTS = [joinpath(@__DIR__, "..", "problems"),
               joinpath(@__DIR__, "..", "test", "CI-runs")]

#
# What each kernel passes, by spatial dimension. Keep in step with:
#   src/kernel/operators/rhs_gpu.jl          (1D _build_rhs_gpu_v0!, 2D/3D BC)
#   src/kernel/operators/rhs_laguerre_gpu.jl (2D Laguerre BC)
#   src/kernel/operators/rhs_jacc.jl         (2D JACC BC)
#
# A case is 1D/2D/3D by how many arguments its own definition takes, so the table
# is keyed on the arity itself: these are the ONLY arities any kernel calls with.
#
const BC_ARITIES = Dict(
    5  => "1D  (q, qe, coords, t, lpert)",
    8  => "2D  (q, qe, coords, t, nx, ny, qbdy, lpert)",
    9  => "3D  (q, qe, coords, t, nx, ny, nz, qbdy, lpert)",
)

"Return (path, arity, paramlist) for every user_bc_dirichlet_gpu definition."
function bc_definitions()
    defs = Tuple{String,Int,String}[]
    for root in ROOTS
        isdir(root) || continue
        for (dir, _, files) in walkdir(root), f in files
            endswith(f, ".jl") || continue
            path = joinpath(dir, f)
            for line in eachline(path)
                m = match(r"^\s*function\s+user_bc_dirichlet_gpu\s*\((.*)\)\s*$", line)
                m === nothing && continue
                params = strip(m.captures[1])
                n = count(==(','), params) + 1
                push!(defs, (relpath(path, joinpath(@__DIR__, "..")), n, params))
            end
        end
    end
    return defs
end

@testset verbose = true "device callback signatures" begin

    defs = bc_definitions()

    @testset "found the cases at all" begin
        # a guard on the guard: if the walk or the regex breaks, the rest of this
        # file would pass vacuously, which is worse than failing
        @test length(defs) > 50
    end

    @testset "user_bc_dirichlet_gpu arity" begin
        bad = filter(d -> !haskey(BC_ARITIES, d[2]), defs)
        if !isempty(bad)
            println("\n  user_bc_dirichlet_gpu definitions the kernels cannot call:")
            for (path, n, params) in bad
                println("    ", path, "\n      ", n, " parameters: (", params, ")")
            end
            println("  accepted arities:")
            for k in sort(collect(keys(BC_ARITIES)))
                println("    ", k, "  ", BC_ARITIES[k])
            end
        end
        @test isempty(bad)
    end

    @testset "the coordinate slot is named coords, not x/y/z" begin
        # The kernels hand over @view(coords[:, ip]) — one cache line per node,
        # which is the reason mesh.coords is laid out (nsd, npoin) in the first
        # place. A parameter still called `x` there is a leftover from before that
        # migration and will silently be a VIEW, not a scalar, if anyone uses it.
        offenders = Tuple{String,String}[]
        for (path, n, params) in defs
            ps = strip.(split(params, ','))
            length(ps) >= 3 || continue
            slot = ps[3]
            slot == "coords" || push!(offenders, (path, slot))
        end
        if !isempty(offenders)
            println("\n  third parameter should be `coords`:")
            for (path, slot) in offenders
                println("    ", path, "  has `", slot, "`")
            end
        end
        @test isempty(offenders)
    end
end
