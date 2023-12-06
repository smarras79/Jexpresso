using LinearSolve
using SnoopCompile
using WriteVTK
import SciMLBase

include("./plotting/jeplots.jl")

#----------------------------------------------------------------------------------------------------------------------------------------------
# ∂q/∂t = RHS -> q(x,t)
#----------------------------------------------------------------------------------------------------------------------------------------------
# PNG
function write_output(SD::NSD_1D, q::Array, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::PNG)

    #Reference values only (definied in initial conditions)
  
    nvar = length(varnames)
    for ivar = 1:nvar
        #plot_results!(SD, mesh, q[1:mesh.npoin,ivar], "initial", OUTPUT_DIR, varnames; iout=1, nvar=nvar, PT=nothing)
        plot_results(SD, mesh, q[1:mesh.npoin,ivar], "initial", OUTPUT_DIR, varnames; iout=1, nvar=nvar, PT=nothing)
    end
end

function write_output(SD::NSD_1D, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::PNG; nvar=1, qexact=zeros(1,nvar), case="")
    
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  "))
   # fig = Figure()
    #colors = ["Blue","Red","Green","Yellow","Black","Purple","Orange"]
    #p = []
    for iout = 1:size(sol.t[:], 1)
        #icolor = mod(iout,size(colors,1))+1
        #color = colors[icolor]
        title = string("sol.u at time ", sol.t[iout])
        #plot_results!(SD, mesh, sol.u[iout][:], title, OUTPUT_DIR, varnames; iout=iout, nvar=nvar, fig=fig,color = color,p=p,PT=nothing)
        plot_results(SD, mesh, sol.u[iout][:], title, OUTPUT_DIR, varnames; iout=iout, nvar=nvar,PT=nothing)
    end
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE ") )
end


function write_output(SD::NSD_1D, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::HDF5; nvar=1, qexact=zeros(1,nvar), case="")
    nothing
end

function write_output(SD::NSD_2D, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::PNG; nvar=1, qexact=zeros(1,nvar), case="")

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

# ASCII
function write_output(SD::NSD_2D, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::ASCII; nvar=1, PT=nothing)
    
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

function write_output(SD::NSD_2D, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::VTK; nvar=1, qexact=zeros(1,nvar), case="")
    
    println(string(" # Writing output to VTK file:", OUTPUT_DIR, "*.vtu ...  ") )
    for iout = 1:size(sol.t[:],1)
        title = @sprintf "Tracer: final solution at t=%6.4f" sol.t[iout]
        write_vtk(SD, mesh, sol.u[iout][:], title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, qexact=qexact, case=case)
    end
    println(string(" # Writing output to VTK file:", OUTPUT_DIR, "*.vtu ... DONE") )
    
end


function write_output(SD::NSD_2D, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::HDF5; nvar=1, qexact=zeros(1,nvar), case="")
    
    println(string(" # Writing restart HDF5 file:", OUTPUT_DIR, "*.h5 ...  ") )
    iout = size(sol.t[:],1)
    title = @sprintf "Final solution at t=%6.4f" sol.t[iout]
    write_hdf5(SD, mesh, sol.u[iout][:], qexact, title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, case=case)
    #end
    println(string(" # Writing restart HDF5 file:", OUTPUT_DIR, "*.h5 ... DONE") )
    
end
function read_output(SD::NSD_2D, INPUT_DIR::String, inputs::Dict, npoin, outformat::HDF5; nvar=1)
    
    println(string(" # Reading restart HDF5 file:", INPUT_DIR, "*.h5 ...  ") )
    q, qe = read_hdf5(SD, INPUT_DIR, inputs, npoin, nvar)
    println(string(" # Reading restart HDF5 file:", INPUT_DIR, "*.h5 ... DONE") )

    return q, qe
end


#------------
# HDF5 writer/reader
#------------
function write_hdf5(SD::NSD_2D, mesh::St_mesh, q::Array, qe::Array, title::String, OUTPUT_DIR::String, inputs::Dict, varnames; iout=1, nvar=1, case="")
    
    #Write one HDF5 file per variable
    for ivar = 1:nvar
        fout_name = string(OUTPUT_DIR, "/var_", ivar, ".h5")
        idx = (ivar - 1)*mesh.npoin
        h5write(fout_name, "q",  q[idx+1:ivar*mesh.npoin]);
        h5write(fout_name, "qe", qe[1:mesh.npoin, ivar]);
    end
end
function read_hdf5(SD::NSD_2D, INPUT_DIR::String, inputs::Dict, npoin, nvar)
    
    q  = zeros(Float64, npoin, nvar+1)
    qe = zeros(Float64, npoin, nvar+1)
    
    #Write one HDF5 file per variable
    for ivar = 1:nvar
        fout_name = string(INPUT_DIR, "/var_", ivar, ".h5")
        idx = (ivar - 1)*npoin
        q[:, ivar]  = convert(Array{Float64, 1}, h5read(fout_name, "q"))
        qe[:, ivar] = convert(Array{Float64, 1}, h5read(fout_name, "qe"))
    end
    
    return q, qe
end
#------------
# VTK writer
#------------
function write_vtk(SD::NSD_2D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String, inputs::Dict, varnames; iout=1, nvar=1, qexact=zeros(1,nvar), case="")
    #nothing
    
    if (mesh.nelem_semi_inf > 0)
        subelem = Array{Int64}(undef, mesh.nelem*(mesh.ngl-1)^2+mesh.nelem_semi_inf*(mesh.ngl-1)*(mesh.ngr-1), 4)
        cells = [MeshCell(VTKCellTypes.VTK_QUAD, [1, 2, 4, 3]) for _ in 1:mesh.nelem*(mesh.ngl-1)^2+mesh.nelem_semi_inf*(mesh.ngl-1)*(mesh.ngr-1)]
    else
        subelem = Array{Int64}(undef, mesh.nelem*(mesh.ngl-1)^2, 4)
        cells = [MeshCell(VTKCellTypes.VTK_QUAD, [1, 2, 4, 3]) for _ in 1:mesh.nelem*(mesh.ngl-1)^2]
    end
    isel = 1
    for iel = 1:mesh.nelem
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
    end
    
    for iel = 1:mesh.nelem_semi_inf
        for i = 1:mesh.ngl-1
            for j = 1:mesh.ngr-1
                ip1 = mesh.connijk_lag[iel,i,j]
                ip2 = mesh.connijk_lag[iel,i+1,j]
                ip3 = mesh.connijk_lag[iel,i+1,j+1]
                ip4 = mesh.connijk_lag[iel,i,j+1]
                subelem[isel, 1] = ip1
                subelem[isel, 2] = ip2
                subelem[isel, 3] = ip3
                subelem[isel, 4] = ip4
                
                cells[isel] = MeshCell(VTKCellTypes.VTK_QUAD, subelem[isel, :])
                
                isel = isel + 1
            end
        end
    end
    
    npoin = mesh.npoin
    qout = copy(q)

    
    if (inputs[:CL] == CL())

        if (inputs[:SOL_VARS_TYPE] == TOTAL())
            
            outvars = ("ρ", "u", "v", "θ")
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
                
                if (size(qexact, 1) === npoin)

                    if inputs[:loutput_pert] == true
                        ("dρ", "u", "v", "dθ")
                        
                        #ρ'
                        qout[1:npoin] .= q[1:npoin] .- qexact[1:npoin,1]
                        
                        #θ' = (ρθ - ρθref)/ρ = ρθ/ρ - ρrefθref/ρref
                        ivar = 4
                        idx = (ivar - 1)*npoin
                        qout[idx+1:4*npoin] .= q[idx+1:4*npoin]./q[1:npoin] .- qexact[1:npoin,4]./qexact[1:npoin,1]
                    else
                        
                        outvars = ("ρ", "u", "v", "θ")
                        ivar = 4
                        idx = (ivar - 1)*npoin
                        qout[idx+1:4*npoin] .= q[idx+1:4*npoin]./q[1:npoin]

                    end
                end
            end
        else
            
            outvars = ("dρ", "u", "v", "dθ")
            
            #ρ
            qout[1:npoin] .= q[1:npoin]
            
            #u = ρu/ρ
            ivar = 2
            idx = (ivar - 1)*npoin
            #qout[idx+1:2*npoin] .= q[idx+1:2*npoin]./(qout[1:npoin] .+ qexact[1:npoin,1])
            qout[idx+1:2*npoin] .= (q[idx+1:2*npoin] .+ qexact[1:npoin,2])./(qout[1:npoin] .+ qexact[1:npoin,1]) .- qexact[1:npoin,2]./qexact[1:npoin,1]
            #v = ρv/ρ
            ivar = 3
            idx = (ivar - 1)*npoin
            #qout[idx+1:3*npoin] .= q[idx+1:3*npoin]./(qout[1:npoin] .+ qexact[1:npoin,1])
            qout[idx+1:3*npoin] .= (q[idx+1:3*npoin] .+ qexact[1:npoin,3])./(qout[1:npoin] .+ qexact[1:npoin,1]) .- qexact[1:npoin,3]./qexact[1:npoin,1]

            if case == "rtb" || case == "mountain"
                
                if (size(qexact, 1) === npoin)
                    
                    ivar = 4
                    idx = (ivar - 1)*npoin
                    qout[idx+1:4*npoin] .= (q[idx+1:4*npoin] .+ qexact[1:npoin,4])./(qout[1:npoin] .+ qexact[1:npoin,1]) .- qexact[1:npoin,4]./qexact[1:npoin,1]
                end
            end
        end
    elseif (inputs[:CL] == NCL())

        outvars = ("ρ", "u", "v", "θ")
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
            
            if (inputs[:loutput_pert] == true && size(qexact, 1) === npoin)
                outvars = ("dρ", "u", "v", "dθ")
                #ρ'
                qout[1:npoin] .= (q[1:npoin] .- qexact[1:npoin,1])
                
                #θ' = (ρθ - ρθref)/ρ = ρθ/ρ - ρrefθref/ρref
                ivar = 4
                idx = (ivar - 1)*npoin
                qout[idx+1:4*npoin] .= (q[idx+1:4*npoin] .- qexact[1:npoin,4])
                
            else
                
                outvars = ("ρ'", "u", "v", "e")
                #E = ρE/ρ
                idx = 4*npoin
                qout[idx+1:4*npoin] .= @views((q[2*npoin+1:4*npoin] .- 0.5*(q[npoin+1:2*npoin].*q[npoin+1:2*npoin] .+ q[npoin+1:3*npoin].*q[npoin+1:3*npoin])./q[1:npoin])./q[1:npoin]) #internal energy: p/((γ-1)ρ)
            end
        end
    end
    #=
    outvars = ("ρ", "ρu", "ρv", "ρθ")
    
    #ρ
    qout[1:npoin] .= q[1:npoin]
    
    #u = ρu/ρ
    ivar = 2
    idx = (ivar - 1)*npoin
    qout[idx+1:2*npoin] .= q[idx+1:2*npoin]

    #v = ρv/ρ
    ivar = 3
    idx = (ivar - 1)*npoin
    qout[idx+1:3*npoin] .= q[idx+1:3*npoin]
    
    ivar = 4
    idx = (ivar - 1)*npoin
    qout[idx+1:4*npoin] .= q[idx+1:4*npoin]
    =#


if nvar > 4
    for ivar = 5:nvar
        idx = (ivar - 1)*npoin
        qout[idx+1:ivar*npoin] .= q[idx+1:ivar*npoin]
    end
end


    #Solution:
    fout_name = string(OUTPUT_DIR, "/iter_", iout, ".vtu")    
    vtkfile = vtk_grid(fout_name, mesh.x[1:npoin], mesh.y[1:npoin], mesh.y[1:npoin]*0.0, cells)
    for ivar = 1:nvar
        idx = (ivar - 1)*npoin
        vtkfile[string(varnames[ivar]), VTKPointData()] =  @view(qout[idx+1:ivar*npoin])
    end
    outfiles = vtk_save(vtkfile)
        
end

function write_vtk_ref(SD::NSD_2D, mesh::St_mesh, q::Array, file_name::String, OUTPUT_DIR::String; iout=1, nvar=1, qexact=zeros(1,nvar), case="", outvarsref=tuple(("" for _ in 1:nvar)))

    #nothing
    if (mesh.nelem_semi_inf > 0)
        subelem = Array{Int64}(undef, mesh.nelem*(mesh.ngl-1)^2+mesh.nelem_semi_inf*(mesh.ngl-1)*(mesh.ngr-1), 4)
        cells = [MeshCell(VTKCellTypes.VTK_QUAD, [1, 2, 4, 3]) for _ in 1:mesh.nelem*(mesh.ngl-1)^2+mesh.nelem_semi_inf*(mesh.ngl-1)*(mesh.ngr-1)]
    else
        subelem = Array{Int64}(undef, mesh.nelem*(mesh.ngl-1)^2, 4)
        cells = [MeshCell(VTKCellTypes.VTK_QUAD, [1, 2, 4, 3]) for _ in 1:mesh.nelem*(mesh.ngl-1)^2]
    end
    
    isel = 1
    for iel = 1:mesh.nelem
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
    end
    
    for iel = 1:mesh.nelem_semi_inf
        for i = 1:mesh.ngl-1
            for j = 1:mesh.ngr-1
                ip1 = mesh.connijk_lag[iel,i,j]
                ip2 = mesh.connijk_lag[iel,i+1,j]
                ip3 = mesh.connijk_lag[iel,i+1,j+1]
                ip4 = mesh.connijk_lag[iel,i,j+1]
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
    
    #Reference values only (definied in initial conditions)
    fout_name = string(OUTPUT_DIR, "/", file_name, ".vtu")
    
    vtkfile = vtk_grid(fout_name, mesh.x[1:mesh.npoin], mesh.y[1:mesh.npoin], mesh.y[1:mesh.npoin]*0.0, cells)

    for ivar = 1:length(outvarsref)
        vtkfile[string(outvarsref[ivar]), VTKPointData()] =  @view(q[1:mesh.npoin,ivar])
    end
    outfiles = vtk_save(vtkfile)
end
