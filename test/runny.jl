using Test
using Jexpresso

function run_example(problem_name::String, case_name::String)
    ENV["JEXPRESSO_HOME"] = joinpath(@__DIR__, "..") 
    example_dir = joinpath(ENV["JEXPRESSO_HOME"], "test","reference", "problems", "equations", problem_name, case_name)
    @testset "$problem_name - $case_name" begin
        cd(example_dir)
        empty!(ARGS) # Clear ARGS to ensure clean state
        push!(ARGS, problem_name, case_name)
        include(joinpath(ENV["JEXPRESSO_HOME"], "src", "Jexpresso.jl"))
    end
end

# Run test sets for each example
@testset "JEXPRESSO Examples" begin
    # List of (problem_name, case_name) tuples
    examples = [
        #("CompEuler", "2d"),  ##errore in q_define sembrebbe initialize ln11
        #("CompEuler", "dc"),   # Initial Mass  :   1.458934183926933e8 and then stops
        ("CompEuler", "dc-mount"), #LoadError: MethodError: no method matching filter!
        #=("CompEuler", "HSmount"),
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
        #("Elliptic", "case1"),=#
    ]

    for (problem_name, case_name) in examples
        run_example(problem_name, case_name)
    end
end
