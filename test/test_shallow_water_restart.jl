# Run from the repository root:
# julia --project=. test/test_shallow_water_restart.jl
# Also supports mpiexec -n 2 with the same command.
using Test, Jexpresso, MPI, PartitionedArrays

const CASE_DIR = joinpath(@__DIR__, "..", "problems", "ShallowWater", "ReducedGravityDSGS")
@eval Jexpresso begin
    ENERGY_EQUATION_THETA = Ref(false)
    parsed_equations = "ShallowWater"
    parsed_equations_case_name = "ReducedGravityDSGS"
    user_input_file = $(joinpath(CASE_DIR, "user_inputs.jl"))
end
for file in ("user_inputs.jl", "initialize.jl", "user_flux.jl", "user_source.jl",
             "user_bc.jl", "user_primitives.jl")
    Base.include(Jexpresso, joinpath(CASE_DIR, file))
end
MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const root = MPI.bcast(rank == 0 ? mktempdir(; prefix="swe_restart_test_") : "", 0, comm)

function deck(directory; restart=false, checkpoint=joinpath(directory, "restart"), tend=2400.0, outformat="hdf5")
    inputs = Jexpresso.user_inputs()
    merge!(inputs, Dict(
        :nop => 2, :tend => tend, :diagnostics_at_times => [0.0, 1200.0, 2400.0],
        :output_dir => directory, :lrestart => restart,
        :restart_time => 1200.0, :restart_input_file_path => checkpoint,
        :restart_output_file_path => joinpath(directory, "restart"),
        :outformat => outformat, :lwrite_initial => false,
        :lprecompile_warmup => false, :lmin_h_diagnostic => false,
        :swe_limiter_diagnostic => false,
    ))
    Jexpresso.mod_inputs_user_inputs!(inputs, rank)
    @eval Jexpresso inputs = $inputs
    mkpath(directory)
    return inputs
end

@testset "ReducedGravityDSGS HDF5 and VTK restart" begin
    with_mpi() do distribute
        continuous = deck(joinpath(root, "continuous"))
        sem, model = Jexpresso.sem_setup(continuous, MPI.Comm_size(comm), distribute)
        function run(inputs)
            q = Jexpresso.initialize(sem.mesh.SD, 0, sem.mesh, inputs, inputs[:output_dir], Float64)
            p, u = Jexpresso.params_setup(sem, q, inputs, inputs[:output_dir], Float64)
            sol = Jexpresso.time_loop!(inputs, p, u, model)
            @test Jexpresso.successful_retcode(sol)
            return sol, p
        end
        first_leg = deck(joinpath(root, "first_leg"); tend=1200.0, outformat="vtk")
        first_leg[:lwrite_initial] = true
        split, before = run(first_leg)
        checkpoint = first_leg[:restart_output_file_path]
        @test !isfile(joinpath(checkpoint, "swe_history_$(rank).h5"))
        # Even an old, unreadable history sidecar must be ignored.
        write(joinpath(checkpoint, "swe_history_$(rank).h5"), "unused legacy history")
        @test before.dsgs_nhist[] >= 3
        @test maximum(abs, split.u[end][sem.mesh.npoin+1:end]) > 0

        resumed = deck(joinpath(root, "resumed"); restart=true, checkpoint)
        restored = Jexpresso.initialize(sem.mesh.SD, 0, sem.mesh, resumed, resumed[:output_dir], Float64)
        @test resumed[:tinit] == 1200.0
        @test vec(restored.qn[:, 1:3]) == split.u[end]
        @test restored.qe[:, 1:3] == before.qp.qe[:, 1:3]
        p, u = Jexpresso.params_setup(sem, restored, resumed, resumed[:output_dir], Float64)
        function fresh_history(p)
            @test p.dsgs_qn == p.qp.qn
            @test p.dsgs_qnm1 == p.qp.qn
            @test p.dsgs_qnm2 == p.qp.qn
            @test p.dsgs_thist[] == -1.0e30
            @test p.dsgs_nhist[] == 0
            @test p.dsgs_wt[] == (0.0, 0.0, 0.0)
            @test !p.dsgs_stage[]
        end
        fresh_history(p)
        restarted = Jexpresso.time_loop!(resumed, p, u, model)
        @test Jexpresso.successful_retcode(restarted)
        @test restarted.t[1] >= 1200.0

        vtk_inputs = deck(joinpath(root, "vtk_resumed"))
        vtk_inputs[:lrestart_vtk] = true
        vtk_inputs[:restart_vtk_input_dir] = first_leg[:output_dir]
        # The snapshot time must not be inferred from this new schedule.
        vtk_inputs[:diagnostics_at_times] = [0.0, 2400.0]
        vtk_q = Jexpresso.initialize(sem.mesh.SD, 0, sem.mesh, vtk_inputs, vtk_inputs[:output_dir], Float64)
        @test vtk_inputs[:tinit] == 1200.0
        @test vtk_inputs[:restart_vtk_iout] == 2
        @test vec(vtk_q.qn[:, 1:3]) == split.u[end]
        @test vtk_q.qe[:, 1:3] == before.qp.qe[:, 1:3]
        vtk_p, vtk_u = Jexpresso.params_setup(sem, vtk_q, vtk_inputs, vtk_inputs[:output_dir], Float64)
        fresh_history(vtk_p)
        vtk_sol = Jexpresso.time_loop!(vtk_inputs, vtk_p, vtk_u, model)
        @test Jexpresso.successful_retcode(vtk_sol)
        err = MPI.Allreduce(maximum(abs.(restarted.u[end] .- vtk_sol.u[end])), MPI.MAX, comm)
        rank == 0 && println("HDF5 versus VTK continuation with fresh history: max difference = ", err)
        @test restarted.u[end] ≈ vtk_sol.u[end] rtol=1e-12 atol=1e-12

        selected = deck(joinpath(root, "vtk_selected"))
        selected[:lrestart_vtk] = true
        selected[:restart_vtk_input_dir] = first_leg[:output_dir]
        selected[:restart_vtk_iout] = 1
        selected[:diagnostics_at_times] = [700.0, 900.0]
        initial = Jexpresso.initialize(sem.mesh.SD, 0, sem.mesh, selected, selected[:output_dir], Float64)
        @test selected[:tinit] == 0.0
        @test initial.qn[:, 1] ≈ fill(150.0, sem.mesh.npoin) atol=1e-12 rtol=0
        @test all(initial.qn[:, 2:3] .== 0.0)

        # Older VTU files have H/Hu/Hv but no dedicated restart fields.
        legacy_dir = joinpath(root, "legacy_vtk")
        mkpath(joinpath(legacy_dir, "iter_2"))
        part = lpad(rank + 1, ndigits(MPI.Comm_size(comm)), '0')
        piece = "iter_2_$(part).vtu"
        bytes = read(joinpath(first_leg[:output_dir], "iter_2", piece))
        appended = findfirst(codeunits("<AppendedData"), bytes)
        gt = findnext(==(UInt8('>')), bytes, last(appended))
        header = replace(String(bytes[1:gt]), "restart_" => "unused_")
        open(joinpath(legacy_dir, "iter_2", piece), "w") do file
            write(file, header)
            write(file, bytes[gt+1:end])
        end
        if rank == 0
            manifest = read(joinpath(first_leg[:output_dir], "iter_2.pvtu"), String)
            write(joinpath(legacy_dir, "iter_2.pvtu"), replace(manifest, "restart_" => "unused_"))
            cp(joinpath(first_leg[:output_dir], "simulation.pvd"), joinpath(legacy_dir, "simulation.pvd"))
        end
        MPI.Barrier(comm)
        legacy = deck(joinpath(root, "legacy_read"))
        legacy[:lrestart_vtk] = true
        legacy[:restart_vtk_input_dir] = legacy_dir
        legacy_q = Jexpresso.initialize(sem.mesh.SD, 0, sem.mesh, legacy, legacy[:output_dir], Float64)
        @test legacy[:tinit] == 1200.0
        @test vec(legacy_q.qn[:, 1:3]) == split.u[end]
        @test all(legacy_q.qe[:, 1] .== 150.0)

        # Checkpoints must also be reached when their cadence is not a
        # multiple of Δt. The last checkpoint here is at 1350, not 900.
        cadence = deck(joinpath(root, "cadence"); tend=1500.0)
        cadence[:restart_time] = 450.0
        cadence[:diagnostics_at_times] = Float64[]
        run(cadence)
        loaded = deck(joinpath(root, "cadence_read"); restart=true,
                      checkpoint=cadence[:restart_output_file_path], tend=1500.0)
        Jexpresso.initialize(sem.mesh.SD, 0, sem.mesh, loaded, loaded[:output_dir], Float64)
        @test loaded[:tinit] == 1350.0

        # The throwaway warm-up extends to 300, but the real solve ends at
        # 150: no checkpoint should be written at all.
        warm = deck(joinpath(root, "warmup"); tend=150.0)
        warm[:restart_time] = 300.0
        warm[:diagnostics_at_times] = Float64[]
        warm[:lprecompile_warmup] = true
        run(warm)
        @test !ispath(warm[:restart_output_file_path])
    end
end
rank == 0 && println("Restart test artifacts: ", root)
