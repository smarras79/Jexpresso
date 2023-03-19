using LinearSolve
using SnoopCompile
using WriteVTK
import SciMLBase

include("./plotting/jeplots.jl")

abstract type AbstractOutFormat end
struct PNG <: AbstractOutFormat end
struct ASCII <: AbstractOutFormat end 

#----------------------------------------------------------------------------------------------------------------------------------------------
# ∂q/∂t = RHS -> q(x,t)
#----------------------------------------------------------------------------------------------------------------------------------------------
# PNG 
function write_output(sol::ODESolution, SD, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG)
    
    for iout = 1: inputs[:ndiagnostics_outputs]
        title = @sprintf "Tracer: final solution at t=%6.4f" sol.t[iout]
        plot_results(SD, mesh.x, mesh.y, sol.u[iout], title, string(OUTPUT_DIR, "/it.", iout, ".png"))
    end
    
end

# ASCII
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



#----------------------------------------------------------------------------------------------------------------------------------------------
# Aq = b -> q(x)
#----------------------------------------------------------------------------------------------------------------------------------------------
#VTK:
function write_vtk(sol::Array, SD, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict)
    cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (i, )) for i = 1:mesh.npoin]
    vtk_grid(string(OUTPUT_DIR, "qsolution"), mesh.x, mesh.y, mesh.z, cells) do vtk
        vtk["q", VTKPointData()] = sol[:,1]
    end
end


# PNG
function write_output(sol::SciMLBase.LinearSolution, SD, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG)
    title = @sprintf "Solution to ∇⋅∇(q) = f"
    plot_results(SD, mesh.x, mesh.y, sol.u, title, string(OUTPUT_DIR, "Axb.png"))
    
end

# ASCII
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_1D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII) nothing end
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_3D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII) nothing end
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_2D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII)    
    title = @sprintf "Solution to ∇⋅∇(q) = f"
    fname = @sprintf "Axb.dat"
    open(string(OUTPUT_DIR, "/", fname), "w") do f
        for ip = 1:length(sol.u)
            @printf(f, " %d %.6f %.6f %.6f \n", ip, mesh.x[ip], mesh.y[ip], sol.u[ip])
        end #f
    end  
end

