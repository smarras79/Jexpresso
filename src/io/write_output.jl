using LinearSolve
using SnoopCompile
using WriteVTK
import SciMLBase

include("./plotting/jeplots.jl")

#----------------------------------------------------------------------------------------------------------------------------------------------
# ∂q/∂t = RHS -> q(x,t)
#----------------------------------------------------------------------------------------------------------------------------------------------
# PNG
function write_output(sol::ODESolution, SD::NSD_1D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG; nvar=1, PT=nothing)
    
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  "))
    for iout = 1:size(sol.t[:], 1)
        title = string("sol.u at time ", sol.t[iout])
        plot_results(SD, mesh, sol.u[iout][:], title, OUTPUT_DIR; iout=iout, nvar=nvar, PT=PT)
    end
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE ") )
end
function write_output(SD::NSD_2D, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG; nvar=1)
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  "))
    
    if inputs[:lplot_surf3d]
        for iout = 1:size(sol.t[:], 1)
            title = @sprintf "Tracer: final solution at t=%6.4f" sol.t[iout]
            plot_surf3d(SD, mesh, sol.u[iout][:], title, OUTPUT_DIR; iout=iout, nvar=nvar, smoothing_factor=inputs[:smoothing_factor])
        end
    else
        for iout = 1:size(sol.t[:],1)
            title = @sprintf "Tracer: final solution at t=%6.4f" sol.t[iout]
            plot_triangulation(SD, mesh, sol.u[iout][:], title,  OUTPUT_DIR; iout=iout, nvar=nvar)
        end
    end
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE"))
end
function write_output(sol::ODESolution, SD::NSD_3D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG; nvar=1)
    nothing
end

# ASCII
function write_output(sol::ODESolution, SD::NSD_1D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII; nvar=1, PT=nothing)

    println(string(" # Writing output to ASCII file:", OUTPUT_DIR, "*.dat ...  ") )
    for iout = 1:size(sol.t[:],1)
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
function write_output(SD::NSD_2D, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII; nvar=1, PT=nothing)
    
    println(string(" # Writing output to ASCII file:", OUTPUT_DIR, "*.dat ...  ") )
    for iout = 1:size(sol.t[:],1)
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

function write_output(SD::NSD_2D, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::VTK; nvar=1, qexact=zeros(1,nvar), case="")
    
    println(string(" # Writing output to VTK file:", OUTPUT_DIR, "*.vtu ...  ") )
    for iout = 1:size(sol.t[:],1)
        title = @sprintf "Tracer: final solution at t=%6.4f" sol.t[iout]
        write_vtk(SD, mesh, sol.u[iout][:], title, OUTPUT_DIR, inputs; iout=iout, nvar=nvar, qexact=qexact, case=case)
    end
    println(string(" # Writing output to VTK file:", OUTPUT_DIR, "*.vtu ... DONE") )
    
end

#----------------------------------------------------------------------------------------------------------------------------------------------
# Aq = b -> q(x)
#----------------------------------------------------------------------------------------------------------------------------------------------
#=
# PNG
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_2D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG; nvar=1)

    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  ") )
    title = @sprintf "Solution to ∇⋅∇(q) = f"
    
    if inputs[:lplot_surf3d]
        plot_surf3d(SD, mesh, sol.u, title, OUTPUT_DIR; iout=1, nvar=1, smoothing_factor=inputs[:smoothing_factor])
    else
        plot_triangulation(SD, mesh, sol.u, title, OUTPUT_DIR;)
    end
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE") )
end

# ASCII
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_1D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict,  outformat::ASCII; nvar=1) nothing end
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_3D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict,  outformat::ASCII; nvar=1) nothing end
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_2D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict,  outformat::ASCII; nvar=1)   
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

function write_output(sol::SciMLBase.LinearSolution, SD::NSD_2D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict,  outformat::VTK; nvar=1)
    
    println(string(" # Writing output to VTK file:", OUTPUT_DIR, "*.vtu ...  ") )
    for iout = 1:inputs[:ndiagnostics_outputs]
        title = @sprintf " ∇²q = f"
        write_vtk(SD, mesh, sol.u, title, OUTPUT_DIR, inputs;)
        #vtkss()
    end
    println(string(" # Writing output to VTK file:", OUTPUT_DIR, "*.vtu ... DONE") )
end
=#
#------------
# VTK writer
#------------
function write_vtk(SD::NSD_2D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String, inputs::Dict; iout=1, nvar=1, qexact=zeros(1,nvar), case="")
    #nothing
    
    subelem = Array{Int64}(undef, mesh.nelem*(mesh.ngl-1)^2, 4)
    cells = [MeshCell(VTKCellTypes.VTK_QUAD, [1, 2, 4, 3]) for _ in 1:mesh.nelem*(mesh.ngl-1)^2]
    
    isel = 1
    for iel = 1:mesh.nelem
        #if iel in mesh.bdy_edge_in_elem[:]
        #    nothing
        #else
            for i = 1:mesh.ngl-1
                for j = 1:mesh.ngl-1
                    ip1 = mesh.connijk[iel,i,j]
                    ip2 = mesh.connijk[iel,i+1,j]
                    ip3 = mesh.connijk[iel,i+1,j+1]
                    ip4 = mesh.connijk[iel,i,j+1]
                    subelem[isel, 1] = ip1
                    subelem[isel, 2] = ip2
                    subelem[isel, 3] = ip3
                    subelem[isel, 4] = ip4
                    
                    cells[isel] = MeshCell(VTKCellTypes.VTK_QUAD, subelem[isel, :])
                    
                    isel = isel + 1
                end
            end
        #end
    end
    
    npoin = mesh.npoin
    qout = copy(q)

    
    if (inputs[:CL] == CL())

        if (inputs[:SOL_VARS_TYPE] == TOTAL())
            #ρ
            qout[1:npoin] .= q[1:npoin]
            
            #u = ρu/ρ
            ivar = 2
            idx = (ivar - 1)*npoin
            qout[idx+1:2*npoin] .= q[idx+1:2*npoin]./q[1:npoin]

            #v = ρv/ρ
            ivar = 3
            idx = (ivar - 1)*npoin
            qout[idx+1:3*npoin] .= q[idx+1:3*npoin]./q[1:npoin]

            if case == "rtb" || case == "mountain"

                outvars = ("rho", "u", "v", "theta")
                
                if (size(qexact, 1) === npoin)

                    if inputs[:loutput_pert] == true
                        
                        #ρ'
                        qout[1:npoin] .= q[1:npoin] .- qexact[1:npoin,1]
                        
                        #θ' = (ρθ - ρθref)/ρ = ρθ/ρ - ρrefθref/ρref
                        ivar = 4
                        idx = (ivar - 1)*npoin
                        qout[idx+1:4*npoin] .= q[idx+1:4*npoin]./q[1:npoin] .- qexact[1:npoin,4]./qexact[1:npoin,1]
                    else
                        ivar = 4
                        idx = (ivar - 1)*npoin
                        qout[idx+1:4*npoin] .= q[idx+1:4*npoin]./q[1:npoin] 
                    end
                end
            end
        else
         
            #ρ
            qout[1:npoin] .= q[1:npoin]
            
            #u = ρu/ρ
            ivar = 2
            idx = (ivar - 1)*npoin
            qout[idx+1:2*npoin] .= q[idx+1:2*npoin]./qout[1:npoin]

            #v = ρv/ρ
            ivar = 3
            idx = (ivar - 1)*npoin
            qout[idx+1:3*npoin] .= q[idx+1:3*npoin]./qout[1:npoin]

            if case == "rtb" || case == "mountain"

                outvars = ("rho", "u", "v", "theta")
                
                if (size(qexact, 1) === npoin)
                    
                    ivar = 4
                    idx = (ivar - 1)*npoin
                    qout[idx+1:4*npoin] .= q[idx+1:4*npoin]./(qout[1:npoin] .+ qexact[1:npoin,1])
                    
                end
            end
        end
    elseif (inputs[:CL] == NCL())
            
        #ρ
        qout[1:npoin] .= (q[1:npoin])
        
        #u = ρu/ρ
        ivar = 2
        idx = (ivar - 1)*npoin
        qout[idx+1:2*npoin] .= (q[idx+1:2*npoin])

        #v = ρv/ρ
        ivar = 3
        idx = (ivar - 1)*npoin
        qout[idx+1:3*npoin] .= (q[idx+1:3*npoin])
        
        if case === "rtb"

            outvars = ("rho", "u", "v", "theta")
            
            if (inputs[:loutput_pert] == true && size(qexact, 1) === npoin)

                #ρ'
                qout[1:npoin] .= (q[1:npoin] .- qexact[1:npoin,1])
                                
                #θ' = (ρθ - ρθref)/ρ = ρθ/ρ - ρrefθref/ρref
                ivar = 4
                idx = (ivar - 1)*npoin
                qout[idx+1:4*npoin] .= (q[idx+1:4*npoin] .- qexact[1:npoin,4])
                                
            else
                outvars = ("rho", "u", "v", "E")
                
                #E = ρE/ρ
                idx = 4*npoin
                qout[idx+1:4*npoin] .= @views((q[2*npoin+1:4*npoin] .- 0.5*(q[npoin+1:2*npoin].*q[npoin+1:2*npoin] .+ q[npoin+1:3*npoin].*q[npoin+1:3*npoin])./q[1:npoin])./q[1:npoin]) #internal energy: p/((γ-1)ρ)
            end
        end
    end
    
    fout_name = string(OUTPUT_DIR, "/rtb_it", iout, ".vtu")
    
    vtkfile = vtk_grid(fout_name, mesh.x[1:npoin], mesh.y[1:npoin], mesh.y[1:npoin]*0.0, cells)
    for ivar = 1:nvar
        idx = (ivar - 1)*npoin
        vtkfile[string(outvars[ivar]), VTKPointData()] =  @view(qout[idx+1:ivar*npoin])
    end
    outfiles = vtk_save(vtkfile)
    
end
