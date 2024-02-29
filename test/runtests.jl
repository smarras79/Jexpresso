using Test
using Jexpresso

function run_example(problem_name::String, case_name::String)
    ENV["JEXPRESSO_HOME"] = joinpath(@__DIR__, "..") 
    example_dir = joinpath(ENV["JEXPRESSO_HOME"], "problems", "equations", problem_name, case_name)
    @testset "$problem_name - $case_name" begin
        cd(example_dir)
        empty!(ARGS) # Clear ARGS to ensure clean state
        push!(ARGS, problem_name, case_name)
        try
            include(joinpath(ENV["JEXPRESSO_HOME"], "src", "Jexpresso.jl"))
            @test true # Passes if no errors occur during execution
        catch e
            @test false # Fails if an error occurs
        end
    end
end

# Run test sets for each example
@testset "JEXPRESSO Examples" begin
    # List of (problem_name, case_name) tuples
    examples = [
        ("CompEuler", "2d"),
        ("CompEuler", "dc"),
        ("CompEuler", "dc-mount"),
        ("CompEuler", "HSmount"),
        ("CompEuler", "HSmount_Lag_working"),
        ("CompEuler", "HSmount_standard"),
        ("CompEuler", "mountain"),
        ("CompEuler", "NHSmount_Lag_working"),
        ("CompEuler", "NHSmount_standard"),
        ("CompEuler", "NLHSmount"),
        ("CompEuler", "nozzleanderson"),
        ("CompEuler", "ScharMount"),
        ("CompEuler", "ScharMount_Lag"),
        ("CompEuler", "theta"),
        ("CompEuler", "theta_laguerre"),
        ("CompEuler", "theta_pert"),
        ("CompEuler", "thetaNC"),
        ("CompEuler", "thetaTracers"),
        ("CompEuler", "wave1d"),
        ("CompEuler", "wave1d_lag"),
        ("AdvDiff", "2d"),
        ("AdvDiff", "2d_Laguerre"),
        ("AdvDiff", "2D_Wave_Train"),
        ("AdvDiff", "case1"),
        ("AdvDiff", "circle"),
        ("AdvDiff", "fd1d"),
        ("AdvDiff", "Simple_Wave"),
        ("AdvDiff", "Wave_Train_Overlapping_Plot"),
        ("AdvDiff", "Wave_Train"),
        ("Helmholtz", "case1"),
        #("Elliptic", "case1"),
    ]

    for (problem_name, case_name) in examples
        run_example(problem_name, case_name)
    end
end
