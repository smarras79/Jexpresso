using Test
using HDF5

# Ensure JEXPRESSO_HOME is set
if !haskey(ENV, "JEXPRESSO_HOME")
    ENV["JEXPRESSO_HOME"] = dirname(dirname(@__DIR__())) # Default path if not set, adjusted to two directories up
end



function find_hdf5_files(directory::String)
    files = readdir(directory, join=true)  # List all files with full paths
    h5_files = filter(f -> occursin(r"\.h5$", f), files)
    if isempty(h5_files)
        error("No HDF5 files found in the directory: $directory")
    end
    return h5_files
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
    # Store the initial directory
    initial_dir = pwd()
    
    try
        ENV["JEXPRESSO_HOME"] = dirname(dirname(@__DIR__()))
        ENV["CI_ENV"] = "true"  # Signal that we are running in the CI environment
        
        example_dir = joinpath(ENV["JEXPRESSO_HOME"],"JExpresso" , "problems", "equations", "CI-runs", parsed_equations, parsed_equations_case_name)
        
        
        test_dir = joinpath(ENV["JEXPRESSO_HOME"],"JExpresso" , "test", "CI-ref", parsed_equations, parsed_equations_case_name)
        test_files = find_hdf5_files(test_dir)

        @testset "$parsed_equations - $parsed_equations_case_name" begin
            cd(example_dir)
            empty!(ARGS)  # Clear ARGS to ensure clean state
            push!(ARGS, parsed_equations, parsed_equations_case_name)
            try
                include(joinpath(ENV["JEXPRESSO_HOME"],"JExpresso" , "src", "Jexpresso.jl"))
            catch e
                error_message = string(e)
                println("Error occurred: ", error_message[1:min(1000, end)])
                @test false
            end
        end

        generated_files = find_hdf5_files(joinpath(ENV["JEXPRESSO_HOME"], "JExpresso" , "problems", "equations", "CI-runs",parsed_equations, parsed_equations_case_name,"output", "CI-runs", parsed_equations, parsed_equations_case_name,"output"))

        for i in 1:length(generated_files)
            generated_file = generated_files[i]
            test_file = test_files[i]  # Assuming the corresponding test file has the same index
    
            test_data = Dict()
            h5open(test_file, "r") do file
                for name in keys(file)
                    test_data[name] = read(file[name])
                end
            end
    
            generated_data = Dict()
            h5open(generated_file, "r") do file
                for name in keys(file)
                    generated_data[name] = read(file[name])
                end
            end 

            @testset "Comparison for reliability check: $(basename(generated_file))" begin
                @test compare_results(test_data, generated_data)
            end
        end

    finally
        # Ensure we navigate back to the initial directory
        cd(initial_dir)
    

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

    for (problem_name, case_name) in examples
        run_example(problem_name, case_name)
    end
end
