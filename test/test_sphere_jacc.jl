#---------------------------------------------------------------------------------
# test/test_sphere_jacc.jl — the portable JACC right-hand side on the CLOSED
#                            SPHERICAL SHELL (src/kernel/operators/sphere_rhs_jacc.jl).
#
#   julia --project=. -t 4 test/test_sphere_jacc.jl
#
# The shell does not go through rhs!, so none of the flat path's tests reach it:
# it has its own residual (sphere_rhs!) and its own time loop. This checks the
# JACC version of that residual against the serial one it replaces, on the grid
# the SWsphere case ships, so the test needs no gmsh binary.
#
# WHAT IS ASSERTED
#
#   * the JACC shell residual equals the serial one to round-off, inviscid and
#     with the Laplace-Beltrami diffusion switched on. Not bit-for-bit: the
#     direct stiffness summation adds the same numbers in a different order (a
#     gather per node instead of a scatter per element), and floating-point
#     addition is not associative;
#   * it is deterministic — two evaluations agree bitwise, whatever the thread
#     count, which is the property the gather buys and an atomic scatter does not;
#   * the residual is not trivially zero, or "they agree" would mean nothing.
#
# The initial state is the Galewsky jet from the case's own initialize.jl, so the
# comparison is made on the field the case actually integrates.
#---------------------------------------------------------------------------------

using Test
using Jexpresso
using Jexpresso: mod_mesh_mesh_driver, build_sphere_metrics, build_sphere_params,
                 build_jacc_cache_sphere, sphere_rhs!
using PartitionedArrays, MPI
using Printf

const CASE = joinpath(@__DIR__, "..", "problems", "ShallowWater", "SWsphere")

# mod_inputs_user_inputs! and the case files read three module globals that
# run.jl normally sets from the command line before including a case.
@eval Jexpresso begin
    parsed_equations           = "ShallowWater"
    parsed_equations_case_name = "SWsphere"
    user_input_file            = $(joinpath(CASE, "user_inputs.jl"))
end
for f in ("user_inputs.jl", "initialize.jl", "user_flux.jl",
          "user_source.jl", "user_bc.jl", "user_primitives.jl")
    Base.include(Jexpresso, joinpath(CASE, f))
end

@testset verbose = true "JACC shell RHS (SWsphere)" begin

    msh = joinpath(CASE, "cubed_sphere.msh")
    if !isfile(msh)
        @warn "skipping: the cubed-sphere grid is not in the repository" msh
    else
        inputs = Jexpresso.user_inputs()
        inputs[:gmsh_filename] = msh
        inputs[:outformat]     = "none"
        Jexpresso.mod_inputs_user_inputs!(inputs, 1)
        inputs[:gmsh_filename] = msh
        @eval Jexpresso global inputs = $inputs

        mesh, metrics, sp, q, SVT = with_mpi() do distribute
            m, _ = mod_mesh_mesh_driver(inputs, 1, distribute)
            mt   = build_sphere_metrics(m, inputs)
            qq   = Jexpresso.initialize(m.SD, nothing, m, inputs, mktempdir(), Jexpresso.TFloat)
            s    = build_sphere_params(m, mt, inputs; neqs = qq.neqs)
            (m, mt, s, qq, inputs[:SOL_VARS_TYPE])
        end

        @info "shell" nelem=Int(mesh.nelem) npoin=Int(mesh.npoin) neqs=sp.neqs lvisc=sp.lvisc threads=Threads.nthreads()

        u   = copy(q.qn)
        neq = sp.neqs

        for lvisc in (false, true)
            @testset "$(lvisc ? "viscous" : "inviscid") matches the serial shell RHS" begin
                μ_save   = copy(sp.μ)
                lv_save  = sp.lvisc
                lvisc || (sp.μ .= 0)
                sp.lvisc = lvisc

                dser = zero(u); djac = zero(u)

                sp.jacc = nothing
                sphere_rhs!(dser, u, q.qe, mesh, metrics, sp, SVT)

                sp.jacc = build_jacc_cache_sphere(mesh, metrics, sp, q.qe, inputs,
                                                  SVT, Jexpresso.TFloat)
                sphere_rhs!(djac, u, q.qe, mesh, metrics, sp, SVT)

                for ieq = 1:neq
                    a = @view(dser[:, ieq]); b = @view(djac[:, ieq])
                    sc = maximum(abs, a)
                    # the Galewsky jet drives every equation; a zero column would
                    # make the comparison below vacuous
                    ieq == 1 || @test sc > 1e-3
                    @test maximum(abs, a .- b) <= 1e-12 * max(sc, eps())
                end

                sp.μ .= μ_save
                sp.lvisc = lv_save
            end
        end

        @testset "bitwise deterministic" begin
            sp.lvisc = true
            r = map(1:2) do _
                d = zero(u)
                sp.jacc = build_jacc_cache_sphere(mesh, metrics, sp, q.qe, inputs,
                                                  SVT, Jexpresso.TFloat)
                sphere_rhs!(d, u, q.qe, mesh, metrics, sp, SVT)
                d
            end
            @test r[1] == r[2]
        end
    end
end
