using Test
using Jexpresso


function run_example(parsed_equations::String, parsed_equations_case_name::String)
    ENV["JEXPRESSO_HOME"] = joinpath(@__DIR__, "..") 
    example_dir = joinpath(ENV["JEXPRESSO_HOME"], "problems", "equations",  parsed_equations, parsed_equations_case_name)
    @testset "$parsed_equations - $parsed_equations_case_name" begin
        cd(example_dir)
        empty!(ARGS) # Clear ARGS to ensure clean state
        push!(ARGS, parsed_equations, parsed_equations_case_name)
        try
            include(joinpath(ENV["JEXPRESSO_HOME"], "src", "Jexpresso.jl"))
            @test true # Passes if no errors occur during execution
        catch e
            error_message = string(e) # Convert error to string
            println("Error occurred: ", error_message[1:min(1000, end)]) # Print out the first xxx characters of the error message(last try 150 needs more to be readable  but it gets overwritten anyway in the github logs)
            @test false # Fails if an error occurs
        end
    end
end

# Run test sets for each example

@testset "JEXPRESSO Examples" begin
    # List of (problem_name, case_name) tuples
    examples = [
        # working
        ("CompEuler", "theta"),
        ("AdvDiff", "2d_Laguerre"),
        ("CompEuler", "dc"),
        ("CompEuler", "nozzleanderson"),
        ("CompEuler", "theta_laguerre"),
        ("CompEuler", "thetaTracers"),
        ("CompEuler", "wave1d"),
        ("CompEuler", "wave1d_lag"),
        ("AdvDiff", "2D_Wave_Train"),
        ("AdvDiff", "case1"),
        ("AdvDiff", "fd1d"),
        ("AdvDiff", "Simple_Wave"),
        ("AdvDiff", "Wave_Train"),
        ("AdvDiff", "Wave_Train_Overlapping_Plot"),
        
        #= not working
        ("AdvDiff", "2d"),
        ("CompEuler", "2d"),
        ("CompEuler", "dc-mount"),
        ("CompEuler", "HSmount"),
        ("CompEuler", "HSmount_Lag_working"),
        ("CompEuler", "HSmount_standard"),
        ("CompEuler", "mountain"),
        ("CompEuler", "NHSmount_Lag_working"),
        ("CompEuler", "NHSmount_standard"),
        ("CompEuler", "NLHSmount"),
        ("CompEuler", "ScharMount"),
        ("CompEuler", "ScharMount_Lag"),      
        ("CompEuler", "theta_pert"),
        ("CompEuler", "thetaNC"),
        
        ("AdvDiff", "circle"),  # OUTDATED     
        ("Elliptic", "case1"),
        ("AdvDiff", "Wave_Train_Overlapping_Plot"),
        ("Helmholtz", "case1"),=#
    ]

    for (problem_name, case_name) in examples
        run_example(problem_name, case_name)
    end
end