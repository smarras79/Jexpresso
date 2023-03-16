include("./plotting/jeplots.jl")

abstract type AbstractOutFormat end
struct PNG <: AbstractOutFormat end
struct ASCII <: AbstractOutFormat end 

function write_output(sol::ODESolution, SD, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG)
    
    for iout = 1: inputs[:ndiagnostics_outputs]
        title = @sprintf "Tracer: final solution at t=%6.4f" sol.t[iout]
        plot_results(SD, mesh.x, mesh.y, sol.u[iout], title, string(OUTPUT_DIR, "/it.", iout, ".png"))
    end
    
end


function write_output(sol::ODESolution, SD::NSD_1D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII)
    
    for iout = 1: inputs[:ndiagnostics_outputs]
        #Write out data at final timestep
	fname = @sprintf "it-%d.dat" iout
    	open(string(OUTPUT_DIR, "/", fname), "w") do f
            for ip = 1:mesh.npoin
                @printf(f, " %d %.6f %.6f \n", ip, mesh.x[ip], sol.u[iout][ip])
            end
        end #f
     
    end
    
end


function write_output(sol::ODESolution, SD::NSD_2D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII)
    
    for iout = 1: inputs[:ndiagnostics_outputs]
        #Write out data at final timestep
	fname = @sprintf "it-%d.dat" iout
    	open(string(OUTPUT_DIR, "/", fname), "w") do f
            for ip = 1:mesh.npoin
                @printf(f, " %d %.6f %.6f %.6f \n", ip, mesh.x[ip], mesh.y[ip], sol.u[iout][ip])
            end
        end #f
     
    end
    
end
