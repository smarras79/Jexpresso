using Test
using Jexpresso

# Define a function to run a test case
function run_example(problem_name::String, case_name::String)
    ENV["JEXPRESSO_HOME"] = joinpath(@__DIR__, "..")  # Update with the actual path
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
        ("CompEuler", "thetaTracers"),
        ("CompEuler", "dc"),
        ("CompEuler", "wave1d"),
        ("CompEuler", "wave1d_lag"),
        ("AdvDiff", "Wave_Train"),
        ("AdvDiff", "Wave_Train_Overlapping_Plot"),
        ("AdvDiff", "2d_Laguerre"),
        ("Helmholtz", "case1"),
        ("CompEuler", "theta_laguerre"),
        ("CompEuler", "HSmount_Lag_working")
    ]

    for (problem_name, case_name) in examples
        run_example(problem_name, case_name)
    end
end
