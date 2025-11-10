using Dates
using DelimitedFiles

final_simulation_time = [0.0]
nop_values = 2:7
flux_names = ["flux_artiano_ec", "flux_artiano_tec", "flux_artiano_etec", "flux_kg", "flux_central"]
T_final = 20.0

# Directory base per i risultati
base_output_dir = "./parametric_runs"
mkpath(base_output_dir)

# File per salvare i risultati delle simulazioni
results_file = joinpath(base_output_dir, "simulation_results.csv")

# Inizializza il file dei risultati
open(results_file, "w") do io
    println(io, "nop,flux_type,final_time,status,timestamp")
end

# Funzione per modificare user_flux.jl
function modify_user_flux!(flux_name::String, flux_file_path::String)
    content = read(flux_file_path, String)
    
    new_function = """
@inline function user_volume_flux(u_ll, u_rr)
\t$flux_name(u_ll, u_rr)
end
"""
    
    pattern = r"@inline function user_volume_flux\(u_ll, u_rr\).*?end"s
    new_content = replace(content, pattern => new_function)
    
    write(flux_file_path, new_content)
    
    println("Modified user_flux.jl to use $flux_name")
end

function modify_user_inputs!(nop_val::Int, output_dir::String, inputs_file_path::String)
    content = read(inputs_file_path, String)
    
    content = replace(content, r":nop\s*=>\s*\d+" => ":nop                 => $nop_val")
    
    content = replace(content, r":output_dir\s*=>\s*\"[^\"]*\"" => ":output_dir          => \"$output_dir\"")
    
    write(inputs_file_path, content)
    
    println("Modified user_inputs.jl: nop=$nop_val, output_dir=$output_dir")
end

problem_dir = "./problems/equations/CompEuler/kelvinHelmholtzChan2022b"
flux_file = joinpath(problem_dir, "user_flux.jl")
inputs_file = joinpath(problem_dir, "user_inputs.jl")

cp(flux_file, flux_file * ".backup", force=true)
cp(inputs_file, inputs_file * ".backup", force=true)
println("Backup files created")

total_runs = length(nop_values) * length(flux_names)

run_counter = 0
for nop_val in nop_values
    for flux_name in flux_names
        global run_counter += 1
        println("\n" * "="^80)
        println("RUN $run_counter/$total_runs: nop=$nop_val, flux=$flux_name")
        println("="^80)
        
        output_dir = joinpath(base_output_dir, "nop$(nop_val)_$(flux_name)")
        mkpath(output_dir)
        
        modify_user_inputs!(nop_val, output_dir, inputs_file)
        modify_user_flux!(flux_name, flux_file)
        
        push!(empty!(ARGS), "CompEuler", "kelvinHelmholtzChan2022b")
        
        timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
        
        include("./src/Jexpresso.jl")
        
	final_time = final_simulation_time[1]
	@show final_time
        status = (final_time >= T_final) ? "completed" : "blow-up"
        
        open(results_file, "a") do io
            println(io, "$nop_val,$flux_name,$final_time,$status,$timestamp")
        end
        
        println("Run completed. Final time: $final_time, Status: $status")
        
        sleep(0.5)
    end
end

# Ripristina i file originali
cp(flux_file * ".backup", flux_file, force=true)
cp(inputs_file * ".backup", inputs_file, force=true)
println("\n" * "="^80)
println("All runs completed! Original files restored.")
println("Results saved in: $results_file")
println("="^80)
