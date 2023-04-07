using LinearSolve
using SnoopCompile
using WriteVTK
import SciMLBase
using Makie

include("./plotting/jeplots.jl")

abstract type AbstractOutFormat end
struct PNG <: AbstractOutFormat end
struct ASCII <: AbstractOutFormat end 

#----------------------------------------------------------------------------------------------------------------------------------------------
# ∂q/∂t = RHS -> q(x,t)
#----------------------------------------------------------------------------------------------------------------------------------------------
# PNG
function write_output(sol::ODESolution, SD::NSD_3D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG; nvar=1) nothing end
function write_output(sol::ODESolution, SD::NSD_2D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG; nvar=1)
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  "))   
    for iout = 1:inputs[:ndiagnostics_outputs]
        title = @sprintf "Tracer: final solution at t=%6.4f" sol.t[iout]
        plot_triangulation(SD, mesh.x, mesh.y, sol.u[iout][:], title,  OUTPUT_DIR; iout=iout, nvar=nvar)
    end
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE"))
end

function write_output(sol::ODESolution, SD::NSD_1D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG; nvar=1)
    
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  "))
    for iout = 1:inputs[:ndiagnostics_outputs]
        title = string("sol.u at time ", sol.t[iout])
        plot_results(SD, mesh.x, mesh.y, sol.u[iout][:], title, OUTPUT_DIR; iout=iout, nvar=nvar)
    end
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.dat ...  DONE ") )
end

# ASCII
function write_output(sol::ODESolution, SD::NSD_1D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII; nvar=1)
    println(string(" # Writing output to ASCII file:", OUTPUT_DIR, "*.dat ...  ") )
    for iout = 1:inputs[:ndiagnostics_outputs]
        title = string("sol.u at time ", sol.t[iout])

        for ivar=1:nvar
            fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".dat")
            open(fout_name, "w") do f
                idx = (ivar - 1)*mesh.npoin
                for ip = (ivar-1)*mesh.npoin+1:ivar*mesh.npoin
                    @printf(f, " %f \n", sol.u[iout][ip])
                end
            end
        end
        #write_ascii(SD, mesh.x, sol.u[iout][:], title, OUTPUT_DIR; iout=iout, nvar=nvar)
    end
    println(string(" # Writing output to ASCII file:", OUTPUT_DIR, "*.dat ...  DONE ") )
end
function write_output(sol::ODESolution, SD::NSD_2D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII; nvar=1)
    
    println(string(" # Writing output to ASCII file:", OUTPUT_DIR, "*.dat ...  ") )
    for iout = 1: inputs[:ndiagnostics_outputs]
        #Write out data at final timestep
	fname = @sprintf "it-%d.dat" iout
    	open(string(OUTPUT_DIR, "/", fname), "w") do f
            for ip = 1:mesh.npoin
                @printf(f, " %d %.6f %.6f %.6f \n", ip, mesh.x[ip], mesh.y[ip], sol.u[iout][ip])
            end
        end #f
    end    
    println(string(" # Writing output to ASCII file:", OUTPUT_DIR, "*.dat ...  DONE ") ) 
end




#----------------------------------------------------------------------------------------------------------------------------------------------
# Aq = b -> q(x)
#----------------------------------------------------------------------------------------------------------------------------------------------
# PNG
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_2D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG; nvar=1)
    
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  ") )
    title = @sprintf "Solution to ∇⋅∇(q) = f"
    plot_triangulation(SD, mesh.x, mesh.y, sol.u, title, OUTPUT_DIR)    
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE") )
    
end

# ASCII
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_1D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII; nvar=1) nothing end
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_3D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII; nvar=1) nothing end
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_2D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII; nvar=1)   
    println(string(" # Writing output to ASCII file:", OUTPUT_DIR, "*.dat ...  ") )
    title = @sprintf "Solution to ∇⋅∇(q) = f"
    fname = @sprintf "Axb.dat"
    open(string(OUTPUT_DIR, "/", fname), "w") do f
        for ip = 1:length(sol.u)
            @printf(f, " %d %.6f %.6f %.6f \n", ip, mesh.x[ip], mesh.y[ip], sol.u[ip])
        end #f
    end
    println(string(" # Writing output to ASCII file:", OUTPUT_DIR, "*.dat ...  DONE") )
end

#VTK:
function write_vtk(sol::Array, SD, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict)
    nothing
    #println(string(" # Writing output to VTK file:", OUTPUT_DIR, "*.vtu ...  ") )
    #cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (i, )) for i = 1:mesh.npoin]
    #vtk_grid(string(OUTPUT_DIR, "qsolution"), mesh.x, mesh.y, mesh.z, cells) do vtk
    #   vtk["q", VTKPointData()] = sol[:,1]
    #end
    #println(string(" # Writing output to VTK file:", OUTPUT_DIR, "*.vtu ...  DONE") )
end
