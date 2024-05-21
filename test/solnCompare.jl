module solnCompare

export run_example

using Test
using HDF5



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
    project_root = dirname(Base.current_project())
    cd(project_root)
    try
        CI_MODE = "true"
        @testset "$parsed_equations - $parsed_equations_case_name" begin
            push!(empty!(ARGS), parsed_equations, parsed_equations_case_name, CI_MODE)
            try
                include(joinpath(project_root , "src", "Jexpresso.jl"))
            catch e
                error_message = string(e)
                println("Error occurred: ", error_message[1:min(1000, end)])
                @test false
            end
        end
        
        ref_dir = joinpath(project_root, "test", "CI-ref", parsed_equations, parsed_equations_case_name, "output")
        output_dir = joinpath(project_root,"test", "CI-runs", parsed_equations, parsed_equations_case_name,"output")
        if !isdir(output_dir) mkpath(output_dir) end
        
        ref_files = find_hdf5_files(ref_dir)
        generated_files = find_hdf5_files(output_dir)

        for i in 1:length(generated_files)
            generated_file = generated_files[i]
            test_file = ref_files[i]  # Assuming the corresponding test file has the same index
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
        nothing
    end
end

end # module
