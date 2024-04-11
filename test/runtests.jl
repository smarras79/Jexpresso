using Test
#using Jexpresso
using HDF5

# Ensure JEXPRESSO_HOME is set
if !haskey(ENV, "JEXPRESSO_HOME")
    ENV["JEXPRESSO_HOME"] = joinpath(@__DIR__, "..")  # Default path if not set
end

function find_first_hdf5_file(directory::String)
    files = readdir(directory, join=true)  # List all files with full paths
    for file in files
        if occursin(r"\.h5$", file)  # Regex to check if the filename ends with .h5
            return file
        end
    end
    error("No HDF5 file found in the directory: $directory")
end

function compare_results(generated_data, expected_data)
    is_equal = true
    for key in keys(generated_data)
        # Compare floating-point numbers up to 8 significant digits
        if typeof(generated_data[key]) == Array{Float64,1}
            is_equal = isapprox(generated_data[key], expected_data[key], atol=1e-8)
        else
            is_equal = generated_data[key] == expected_data[key]
        end
        if !is_equal
            break
        end
    end
    return is_equal
end

function run_example(parsed_equations::String, parsed_equations_case_name::String)
    example_dir = joinpath(ENV["JEXPRESSO_HOME"],  "problems","equations", "CI-runs", parsed_equations, parsed_equations_case_name)
    expected_output_dir = joinpath(ENV["JEXPRESSO_HOME"], "test", "CI-ref", parsed_equations, parsed_equations_case_name)
    expected_output_file = find_first_hdf5_file(expected_output_dir)
    
    expected_data = Dict()
    h5open(expected_output_file, "r") do file
        for name in keys(file)
            expected_data[name] = read(file[name])
        end
    end

    @testset "$parsed_equations - $parsed_equations_case_name" begin
        cd(example_dir)
        empty!(ARGS)  # Clear ARGS to ensure clean state
        push!(ARGS, parsed_equations, parsed_equations_case_name)
        try
            include(joinpath(ENV["JEXPRESSO_HOME"], "src", "Jexpresso.jl"))
            generated_file = find_first_hdf5_file(joinpath(ENV["JEXPRESSO_HOME"], "output","CompEuler","theta","output"))

            generated_data = Dict()
            h5open(generated_file, "r") do file
                for name in keys(file)
                    generated_data[name] = read(file[name])
                end
            end

            @test compare_results(generated_data, expected_data)
        catch e
            error_message = string(e)
            println("Error occurred: ", error_message[1:min(1000, end)])
            @test false
        end
    end
end

# Run test sets for each example
@testset "JEXPRESSO Examples" begin
    examples = [
        # working
        ("CompEuler", "theta"),
        #=("AdvDiff", "2d_Laguerre"),
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
        =#
        
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
    results_dir = ENV["JEXPRESSO_HOME"]

    for (problem_name, case_name) in examples
        run_example(problem_name, case_name)
    end
end
