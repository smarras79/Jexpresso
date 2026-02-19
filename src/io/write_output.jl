using WriteVTK
using Base64

include("./plotting/jeplots.jl")

#------------------------------------------------------------------
# Callback for missing user_uout!()
#------------------------------------------------------------------
function call_user_uout(uout, u, qe, mp, ET, npoin, nvar, noutvar)
    
    if function_exists(@__MODULE__, :user_uout!)
        for ip=1:npoin
            user_uout!(ip, ET, @view(uout[ip,1:noutvar]), @view(u[ip,1:nvar]), @view(qe[ip,1:nvar]); mp=mp)
        end
    else
        for ip=1:npoin
            callback_user_uout!(@view(uout[ip,1:noutvar]), @view(u[ip,1:nvar]))
        end
    end
end

@inline function callback_user_uout!(uout, usol)
    uout[1:end] = usol[1:end]
end

function function_exists(module_name::Module, function_name::Symbol)
    return isdefined(module_name, function_name) && isa(getfield(module_name, function_name), Function)
end
#------------------------------------------------------------------
# END Callback for missing user_uout!()
#------------------------------------------------------------------

function write_output(SD::NSD_1D, q::Array, t, iout, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::PNG; nvar=1, qexact=zeros(1,nvar), case="")
    #OK
    nvar = length(varnames)
    qout = zeros(mesh.npoin)
    
    plot_results(SD, mesh, q[:], "initial", OUTPUT_DIR, varnames, inputs; iout=1, nvar=nvar, PT=nothing)
end

function write_output(SD::NSD_1D, sol, uaux, t, iout,  mesh::St_mesh, mp, 
                      connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                      OUTPUT_DIR::String, inputs::Dict,
                      varnames, outvarnames,
                      outformat::PNG;
                      nvar=1, qexact=zeros(1,nvar), case="")
        
    #
    # 1D PNG of q(t) from dq/dt = RHS
    #
    if (inputs[:plot_overlap])
        fig = nothing  # Plot figure will be created in plot_results!
        colors = ["Blue","Red","Green","Yellow","Black","Purple","Orange"]
        markers = [:circle, :rect, :diamond,:hexagon,:cross,:xcross,:utriangle,:dtriangle,:pentagon,:star4,:star8]
        p = []
        #for iout = 1:size(sol.t[:], 1)
            icolor = mod(iout,size(colors,1))+1
            color = colors[icolor]
            imarker = mod(iout,size(markers,1))+1
            marker = markers[imarker]
            title = string("sol.u at time ", t)
            if (inputs[:backend] == CPU())
                plot_results!(SD, mesh, sol, title, OUTPUT_DIR, varnames, inputs; iout=iout, nvar=nvar, fig=fig,color = color,p=p,marker=marker,PT=nothing)
            else
                uout = KernelAbstractions.allocate(CPU(),Float32, Int64(mesh.npoin))
                KernelAbstractions.copyto!(CPU(), uout, sol)
                plot_results!(SD, mesh, uout, title, OUTPUT_DIR, varnames, inputs; iout=iout, nvar=nvar, fig=fig,color = color,p=p,marker=marker,PT=nothing)
            end
        #end
    else
        #for iout = 1:size(sol.t[:], 1)
        title = string("sol at time ", t)
            if (inputs[:backend] == CPU())
                plot_results(SD, mesh, sol, title, OUTPUT_DIR, varnames, inputs; iout=iout, nvar=nvar,PT=nothing)
            else
                uout = KernelAbstractions.allocate(CPU(), TFloat, Int64(mesh.npoin*nvar))
                KernelAbstractions.copyto!(CPU(), uout, sol)
                convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
                plot_results(SD, mesh, uout, title, OUTPUT_DIR, varnames, inputs; iout=iout, nvar=nvar,PT=nothing)
            end
        #end
    end
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE ") )
end


function write_output(SD, sol::SciMLBase.LinearSolution, uaux, mesh::St_mesh,
                      OUTPUT_DIR::String, inputs::Dict,
                      varnames, outvarnames,
                      outformat::VTK;
                      nvar=1, qexact=zeros(1,nvar), case="")

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)

    #
    # 2D VTK of x from Ax=b
    #    
    if (inputs[:backend] == CPU())

        title = @sprintf "Solution-Axb"
        write_vtk(SD, mesh, sol.u, uaux, nothing, 
                  nothing, nothing,
                  0.0, 0.0, 0.0, 0.0, title, OUTPUT_DIR, inputs,
                  varnames, outvarnames;
                  iout=1, nvar=nvar, qexact=qexact, case=case) 
        
    else
        u = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin*nvar)
        KernelAbstractions.copyto!(CPU(),u,sol.u)

        u_exact = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin,nvar+1)
        KernelAbstractions.copyto!(CPU(),u_exact,qexact)

        convert_mesh_arrays_to_cpu!(SD, mesh, inputs)        
        title = @sprintf "Solution"
        write_vtk(SD, mesh, u, "1", title, OUTPUT_DIR, inputs, varnames; iout=1, nvar=nvar, qexact=u_exact, case=case)
    end
    
    println(string(" # Writing output to VTK file:", OUTPUT_DIR, "*.vtu ... DONE") )
    
end


function write_output(SD, sol, uaux, t, iout,  mesh::St_mesh, mp, 
                      connijk_original, poin_in_bdy_face_original,
                      x_original, y_original, z_original,
                      OUTPUT_DIR::String, inputs::Dict,
                      varnames, outvarnames,
                      outformat::VTK;
                      nvar=1, qexact=zeros(1,nvar), case="")
    
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    title = @sprintf "final solution at t=%6.4f" iout
    if (inputs[:backend] == CPU())

        write_vtk(SD, mesh, sol, uaux, mp, 
                  connijk_original, poin_in_bdy_face_original,
                  x_original, y_original, z_original,
                  t, title, OUTPUT_DIR, inputs,
                  varnames, outvarnames;
                  iout=iout, nvar=nvar, qexact=qexact, case=case) 
        
    else
        #VERIFY THIS on GPU
        u = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin*(nvar+1))
        KernelAbstractions.copyto!(CPU(),u, u_sol)
        u_exact = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin,nvar+1)
        KernelAbstractions.copyto!(CPU(),u_exact,qexact)
        convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
        write_vtk(SD, mesh, u, mp, t, title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, qexact=u_exact, case=case)
    end

    println_rank(string(" # writing ", OUTPUT_DIR, "/iter", iout, ".vtu at t=", t, " s... DONE"); msg_rank = rank )

end

function write_output(SD, sol, uaux, t, iout,  mesh::St_mesh, mp, 
                    connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                    OUTPUT_DIR::String, inputs::Dict,
                    varnames, outvarnames,
                    outformat::NETCDF;
                    nvar=1, qexact=zeros(1,nvar), case="")

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    title = @sprintf "final solution at t=%6.4f" iout
    if (inputs[:backend] == CPU())

        write_NetCDF(SD, mesh, sol, uaux, mp, 
                     connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                     t, title, OUTPUT_DIR, inputs,
                     varnames, outvarnames;
                     iout=iout, nvar=nvar, qexact=qexact, case=case) 

    else
        #VERIFY THIS on GPU
        u = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin*(nvar+1))
        KernelAbstractions.copyto!(CPU(),u, u_sol)
        u_exact = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin,nvar+1)
        KernelAbstractions.copyto!(CPU(),u_exact,qexact)
        convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
        write_NetCDF(SD, mesh, u, uaux, mp, 
                     connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                     t, title, OUTPUT_DIR, inputs,
                     varnames, outvarnames;
                     iout=iout, nvar=nvar, qexact=u_exact, case=case)
    end

    println_rank(string(" # writing ", OUTPUT_DIR, "/iter", iout, ".vtu at t=", t, " s... DONE"); msg_rank = rank )

end

#------------
# VTK writer
#------------
function write_vtk(SD::NSD_2D, mesh::St_mesh, q::Array, qaux::Array, mp, 
                   connijk_original, poin_in_bdy_face_original,
                   x_original, y_original, z_original,
                   t, title::String, OUTPUT_DIR::String, inputs::Dict,
                   varnames, outvarnames;
                   iout=1, nvar=1, qexact=zeros(1,nvar), case="")

    if (isa(varnames, Tuple)    || isa(varnames, String) )   varnames    = collect(varnames) end
    if (isa(outvarnames, Tuple) || isa(outvarnames, String)) outvarnames = collect(outvarnames) end
    
    nvar     = size(varnames, 1)
    noutvar  = size(outvarnames,1) #max(nvar, size(outvarnames,1))
    new_size = size(mesh.x,1)
    if (mesh.nelem_semi_inf > 0)
        subelem = Array{Int64}(undef, mesh.nelem*(mesh.ngl-1)^2+mesh.nelem_semi_inf*(mesh.ngl-1)*(mesh.ngr-1), 4)
        cells = [MeshCell(VTKCellTypes.VTK_QUAD, [1, 2, 4, 3]) for _ in 1:mesh.nelem*(mesh.ngl-1)^2+mesh.nelem_semi_inf*(mesh.ngl-1)*(mesh.ngr-1)]
    else
        subelem = Array{Int64}(undef, mesh.nelem*(mesh.ngl-1)^2, 4)
        cells = [MeshCell(VTKCellTypes.VTK_QUAD, [1, 2, 4, 3]) for _ in 1:mesh.nelem*(mesh.ngl-1)^2]
    end
    isel = 1
    npoin = mesh.npoin
    conn = zeros(mesh.nelem,mesh.ngl,mesh.ngl)
    conn .= mesh.connijk
    
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

    #
    # Fetch user-defined diagnostic vars or take them from the solution vars:
    #
    qout = zeros(Float64, npoin, noutvar)
    u2uaux!(qaux, q, nvar, npoin)
    call_user_uout(qout, qaux, qexact, mp, inputs[:SOL_VARS_TYPE], npoin, nvar, noutvar)

    
    #
    # Write solution to vtk:
    #
    fout_name = string(OUTPUT_DIR, "/iter_", iout)
    vtkfile = map(mesh.parts) do part
        vtkf = pvtk_grid(fout_name,
                         mesh.x[1:mesh.npoin],
                         mesh.y[1:mesh.npoin],
                         mesh.y[1:mesh.npoin]*TFloat(0.0),
                         cells,
                         compress=false;
                         part=part, nparts=mesh.nparts, ismain=(part==1))
        vtkf["part", VTKCellData()] = ones(isel -1) * part
        # vtkf["gid", VTKCellData()] = mesh.ip2gip[:]

        for ivar = 1:noutvar
            idx = (ivar - 1)*npoin
            vtkf[string(outvarnames[ivar]), VTKPointData()] = @view(qout[1:npoin,ivar])
        end
        
        vtkf
    end
    
    outfiles = map(vtk_save, vtkfile)
    
end

function write_vtk(SD::NSD_3D, mesh::St_mesh, q::Array, qaux::Array, mp, 
                   connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                   t, title::String, OUTPUT_DIR::String, inputs::Dict, varnames, outvarnames;
                   iout=1, nvar=1, qexact=zeros(1,nvar), case="")

    if (isa(varnames, Tuple)    || isa(varnames, String) )   varnames    = collect(varnames) end
    if (isa(outvarnames, Tuple) || isa(outvarnames, String)) outvarnames = collect(outvarnames) end
    
    nvar    = size(varnames, 1)
    noutvar = size(outvarnames,1) #max(nvar, size(outvarnames,1))
    npoin   = mesh.npoin
    
    xx = zeros(size(mesh.x,1))
    yy = zeros(size(mesh.x,1))
    zz = zeros(size(mesh.x,1))
    xx .= mesh.x 
    yy .= mesh.y
    zz .= mesh.z
    conn = zeros(mesh.nelem,mesh.ngl,mesh.ngl,mesh.ngl)
    conn .= mesh.connijk
    x_spare = zeros(Bool,size(mesh.x,1),1)
    y_spare = zeros(Bool,size(mesh.y,1),1)
    z_spare = zeros(Bool,size(mesh.z,1),1)
    connijk_spare = zeros(mesh.nelem,mesh.ngl,mesh.ngl,mesh.ngl)
    poin_bdy = zeros(size(mesh.bdy_face_type,1),mesh.ngl,mesh.ngl)
    poin_bdy .= mesh.poin_in_bdy_face
   
    subelem = Array{Int64}(undef, mesh.nelem*(mesh.ngl-1)^3, 8)
    cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, [1, 2, 3, 4, 5, 6, 7, 8]) for _ in 1:mesh.nelem*(mesh.ngl-1)^3]
    
    isel = 1
    for iel = 1:mesh.nelem
        for i = 1:mesh.ngl-1
            for j = 1:mesh.ngl-1
                for k = 1:mesh.ngl-1
                    ip1 = mesh.connijk[iel,i,j,k]
                    ip2 = mesh.connijk[iel,i+1,j,k]
                    ip3 = mesh.connijk[iel,i+1,j+1,k]
                    ip4 = mesh.connijk[iel,i,j+1,k]
                    
                    ip5 = mesh.connijk[iel,i,j,k+1]
                    ip6 = mesh.connijk[iel,i+1,j,k+1]
                    ip7 = mesh.connijk[iel,i+1,j+1,k+1]
                    ip8 = mesh.connijk[iel,i,j+1,k+1]

                    subelem[isel, 1] = ip1
                    subelem[isel, 2] = ip2
                    subelem[isel, 3] = ip3
                    subelem[isel, 4] = ip4
                    subelem[isel, 5] = ip5
                    subelem[isel, 6] = ip6
                    subelem[isel, 7] = ip7
                    subelem[isel, 8] = ip8
                    
                    cells[isel] = MeshCell(VTKCellTypes.VTK_HEXAHEDRON, subelem[isel, :])
                    
                    isel = isel + 1
                end
            end
        end
    end
    
    #
    # Fetch user-defined diagnostic vars or take them from the solution vars:
    #
    qout = zeros(Float64, npoin, noutvar)
    u2uaux!(qaux, q, nvar, npoin)
    call_user_uout(qout, qaux, qexact, mp, inputs[:SOL_VARS_TYPE], npoin, nvar, noutvar)
    
    
    #
    # Write solution:
    #
    fout_name = string(OUTPUT_DIR, "/iter_", iout)
    vtkfile = map(mesh.parts) do part
        vtkf = pvtk_grid(fout_name,
                         mesh.x[1:mesh.npoin],
                         mesh.y[1:mesh.npoin],
                         mesh.z[1:mesh.npoin],
                         cells,
                         compress=false;
                         part=part, nparts=mesh.nparts, ismain=(part==1))
        vtkf["part", VTKCellData()] = ones(isel -1) * part

        for ivar = 1:noutvar
            idx = (ivar - 1)*npoin
            vtkf[string(outvarnames[ivar]), VTKPointData()] = @view(qout[1:npoin,ivar])
        end
        
        vtkf
    end
    
    outfiles = map(vtk_save, vtkfile)
    
end

function write_vtk_grid_only(SD::NSD_2D, mesh::St_mesh, file_name::String, OUTPUT_DIR::String, parts, nparts)
    
    #nothing
    subelem = Array{Int64}(undef, mesh.nelem*(mesh.ngl-1)^2, 4)
    cells = [MeshCell(VTKCellTypes.VTK_QUAD, [1, 2, 4, 3]) for _ in 1:mesh.nelem*(mesh.ngl-1)^2]
    
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
    
    vtkfile = map(parts) do part
        vtkf = pvtk_grid(file_name, mesh.x[1:mesh.npoin], mesh.y[1:mesh.npoin], mesh.y[1:mesh.npoin]*TFloat(0.0), cells, compress=false;
                        part=part, nparts=nparts, ismain=(part==1))
        vtkf["part", VTKCellData()] = ones(isel -1) * part
        vtkf
    end

    
    outfiles = map(vtk_save, vtkfile)
end



function write_vtk_grid_only(SD::NSD_3D, mesh::St_mesh, file_name::String, OUTPUT_DIR::String, parts, nparts)

    subelem = Array{Int64}(undef, mesh.nelem*(mesh.ngl-1)^3, 8)
    cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, [1, 2, 3, 4, 5, 6, 7, 8]) for _ in 1:mesh.nelem*(mesh.ngl-1)^3]
        
    isel = 1
    for iel = 1:mesh.nelem
        for i = 1:mesh.ngl-1
            for j = 1:mesh.ngl-1
                for k = 1:mesh.ngl-1
                    ip1 = mesh.connijk[iel,i,j,k]
                    ip2 = mesh.connijk[iel,i+1,j,k]
                    ip3 = mesh.connijk[iel,i+1,j+1,k]
                    ip4 = mesh.connijk[iel,i,j+1,k]
                    
                    ip5 = mesh.connijk[iel,i,j,k+1]
                    ip6 = mesh.connijk[iel,i+1,j,k+1]
                    ip7 = mesh.connijk[iel,i+1,j+1,k+1]
                    ip8 = mesh.connijk[iel,i,j+1,k+1]

                    subelem[isel, 1] = ip1
                    subelem[isel, 2] = ip2
                    subelem[isel, 3] = ip3
                    subelem[isel, 4] = ip4
                    subelem[isel, 5] = ip5
                    subelem[isel, 6] = ip6
                    subelem[isel, 7] = ip7
                    subelem[isel, 8] = ip8
                    
                    cells[isel] = MeshCell(VTKCellTypes.VTK_HEXAHEDRON, subelem[isel, :])
                    
                    isel = isel + 1
                end
            end
        end
    end
    
    #Reference values only (definied in initial conditions)
    fout_name = string(OUTPUT_DIR, "/", file_name, ".vtu")
    
    # vtkfile = vtk_grid(fout_name, mesh.x[1:mesh.npoin], mesh.y[1:mesh.npoin], mesh.y[1:mesh.npoin]*TFloat(0.0), cells)
    vtkfile = map(parts) do part
        vtkf = pvtk_grid(file_name, mesh.x[1:mesh.npoin], mesh.y[1:mesh.npoin], mesh.z[1:mesh.npoin], cells, compress=false;
                        part=part, nparts=nparts, ismain=(part==1))
        vtkf["part", VTKCellData()] = ones(isel -1) * part
        vtkf
    end
    outfiles = map(vtk_save, vtkfile)
    # outfiles = vtk_save(vtkfile)
end

#------------
# HDF5 writer/reader
#------------
function write_output(SD, sol, uaux, t, iout,  mesh::St_mesh, mp, 
                      connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                      OUTPUT_DIR::String, inputs::Dict,
                      varnames, outvarnames,
                      outformat::HDF5;
                      nvar=1, qexact=zeros(1,nvar), case="")
    
    # println(string(" # Writing restart HDF5 file:", OUTPUT_DIR, "*.h5 ...  ") )
    iout = size(t,1)
    title = @sprintf "Final solution at t=%6.4f" t
    if !isdir(OUTPUT_DIR)
        mkpath(OUTPUT_DIR)
    end
    if (inputs[:backend] == CPU())
    
        write_hdf5(SD, mesh, sol, qexact, t, title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, case=case)
    else
        u_gpu = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin*nvar)
        KernelAbstractions.copyto!(CPU(),u_gpu, u)
        u_exact = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin,nvar+1)
        KernelAbstractions.copyto!(CPU(),u_exact,qexact)
        convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
        write_hdf5(SD, mesh, u_gpu, u_exact, title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, case=case)
        convert_mesh_arrays!(SD, mesh, inputs[:backend], inputs)
    end

    
    # println(string(" # Writing restart HDF5 file:", OUTPUT_DIR, "*.h5 ... DONE") )
    
end
function write_output(SD, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::HDF5; nvar=1, qexact=zeros(1,nvar), case="")
    
    #println(string(" # Writing restart HDF5 file:", OUTPUT_DIR, "*.h5 ...  ") )
    
    iout = size(sol.t[:],1)
    title = @sprintf "Final solution at t=%6.4f" sol.t[iout]

    if (inputs[:backend] == CPU())
        write_hdf5(SD, mesh, sol.u[iout][:], qexact, title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, case=case)
    else
        u_gpu = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin*nvar)
        KernelAbstractions.copyto!(CPU(),u_gpu, sol.u[iout][:])
        u_exact = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin,nvar+1)
        KernelAbstractions.copyto!(CPU(),u_exact,qexact)
        convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
        write_hdf5(SD, mesh, u_gpu, u_exact, title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, case=case)
        convert_mesh_arrays!(SD, mesh, inputs[:backend], inputs)
    end
    
    println(string(" # Writing restart HDF5 file:", OUTPUT_DIR, "*.h5 ... DONE") )
    
end
function read_output(SD, INPUT_DIR::String, inputs::Dict, npoin, outformat::HDF5; nvar=1)
    
    #println(string(" # Reading restart HDF5 file:", INPUT_DIR, "*.h5 ...  ") )
    q, qe = read_hdf5(SD, INPUT_DIR, inputs, npoin, nvar)
    println_rank(string(" # Reading restart HDF5 file:", INPUT_DIR, "*.h5 ... DONE") ; msg_rank = rank)

    return q, qe
end


function write_hdf5(SD, mesh::St_mesh, q::AbstractArray, qe::AbstractArray, t, title::String, OUTPUT_DIR::String, inputs::Dict, varnames; iout=1, nvar=1, case="")
    
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    mpi_size = MPI.Comm_size(comm)
    #Write one HDF5 file timestep
    if rank == 0
        fout_name = string(OUTPUT_DIR, "/t.h5")
        h5open(fout_name, "w") do fid        
            write(fid, "time",  t);
        end
    end
    #Write one HDF5 file per variable
    for ivar = 1:nvar
        fout_name = string(OUTPUT_DIR, "/var_", ivar,"_",rank, ".h5")
        idx = (ivar - 1)*mesh.npoin
        
        h5open(fout_name, "w") do fid        
            write(fid, "q",  q[idx+1:ivar*mesh.npoin]);
            write(fid, "qe", qe[1:mesh.npoin, ivar]);
        end

    end
end

function read_hdf5(SD, INPUT_DIR::String, inputs::Dict, npoin, nvar)
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    mpi_size = MPI.Comm_size(comm)

    q  = zeros(Float64, npoin, nvar+1)
    qe = zeros(Float64, npoin, nvar+1)
    

    #read one HDF5 file time
    
    fout_name = string(INPUT_DIR, "/t.h5")
    time = rank == 0 ? convert(Float64, h5read(fout_name, "time")) : 0.0
    time = MPI.bcast(time, 0, comm)
    inputs[:tinit] = time
    #Write one HDF5 file per variable
    for ivar = 1:nvar
        fout_name = string(INPUT_DIR, "/var_", ivar,"_",rank, ".h5")
        idx = (ivar - 1)*npoin
        q[:, ivar]  = convert(Array{Float64, 1}, h5read(fout_name, "q"))
        qe[:, ivar] = convert(Array{Float64, 1}, h5read(fout_name, "qe"))
    end
    
    return q, qe
end

function write_NetCDF(SD::NSD_2D, mesh::St_mesh, q::Array, qaux::Array, mp, 
                   connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                   t, title::String, OUTPUT_DIR::String, inputs::Dict, varnames, outvarnames;
                   iout=1, nvar=1, qexact=zeros(1,nvar), case="")

    if (isa(varnames, Tuple)    || isa(varnames, String) )   varnames    = collect(varnames) end
    if (isa(outvarnames, Tuple) || isa(outvarnames, String)) outvarnames = collect(outvarnames) end
    
    xx = zeros(size(mesh.x,1))
    yy = zeros(size(mesh.x,1))
    xx .= mesh.x 
    yy .= mesh.y
    nvar     = size(varnames, 1)
    noutvar  = max(nvar, size(outvarnames,1))
    if (mesh.nelem_semi_inf > 0)
        subelem = Array{Int64}(undef, mesh.nelem*(mesh.ngl-1)^2+mesh.nelem_semi_inf*(mesh.ngl-1)*(mesh.ngr-1), 4)
    else
        subelem = Array{Int64}(undef, mesh.nelem*(mesh.ngl-1)^2, 4)
    end
    nsubelem = size(subelem, 1)
    isel = 1
    npoin = mesh.npoin
    conn = zeros(mesh.nelem,mesh.ngl,mesh.ngl)
    conn .= mesh.connijk
    
    for iel = 1:mesh.nelem
        for i = 1:mesh.ngl-1
            for j = 1:mesh.ngl-1
                ip1 = mesh.connijk[iel,i,j]
                ip2 = mesh.connijk[iel,i+1,j]
                ip3 = mesh.connijk[iel,i+1,j+1]
                ip4 = mesh.connijk[iel,i,j+1]
                subelem[isel, 1] = mesh.ip2gip[ip1]
                subelem[isel, 2] = mesh.ip2gip[ip2]
                subelem[isel, 3] = mesh.ip2gip[ip3]
                subelem[isel, 4] = mesh.ip2gip[ip4]
                
                
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
                subelem[isel, 1] = mesh.ip2gip[ip1]
                subelem[isel, 2] = mesh.ip2gip[ip2]
                subelem[isel, 3] = mesh.ip2gip[ip3]
                subelem[isel, 4] = mesh.ip2gip[ip4]
                
                
                isel = isel + 1
            end
        end
    end

    #
    # Fetch user-defined diagnostic vars or take them from the solution vars:
    #
    qout = zeros(Float64, npoin, noutvar)
    u2uaux!(qaux, q, nvar, npoin)
    call_user_uout(qout, qaux, qexact, mp, inputs[:SOL_VARS_TYPE], npoin, nvar, noutvar)


    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    
    local_list = findall(x->x == rank, mesh.gip2owner)
    # Gather data from all processes
    all_ip2gip  = MPI.gather(mesh.ip2gip[local_list], comm)
    all_xx      = MPI.gather(xx[local_list], comm)
    all_yy      = MPI.gather(yy[local_list], comm)
    all_subelem = MPI.gather(subelem, comm)
    all_qout    = MPI.gather(qout[local_list,:], comm)
    
    # Only rank 0 writes the file
    if rank == 0

        # Reshape all data
        global_ip2gip       = vcat(all_ip2gip...)
        global_xx           = vcat(all_xx...)[global_ip2gip]
        global_yy           = vcat(all_yy...)[global_ip2gip]
        global_qout_flat    = vcat(all_qout...)
        global_subelem_flat = vcat(all_subelem...)
        
        global_npoin::Int64    = length(global_xx)
        global_nsubelem::Int64 = length(global_subelem_flat) / 4

        global_subelem = reshape(global_subelem_flat, global_nsubelem, 4)
        global_qout    = reshape(global_qout_flat, global_npoin, noutvar)[global_ip2gip,:]
        
        # Write NetCDF file
        fout_name_global = string(OUTPUT_DIR, "/iter_", iout, ".nc")
        NCDataset(fout_name_global, "c") do ds
            
            # dimensions
            defDim(ds, "nMesh2_node", global_npoin)
            defDim(ds, "nMesh2_face", global_nsubelem)
            defDim(ds, "nMaxMesh2_face_nodes", 4)

            # mesh topology
            mesh = defVar(ds, "mesh", Int32, ())
            mesh.attrib["cf_role"]                = "mesh_topology"
            mesh.attrib["long_name"]              = "Topology of a 2-d unstructured mesh"
            mesh.attrib["topology_dimension"]     = 2
            mesh.attrib["node_coordinates"]       = "Mesh2_node_x Mesh2_node_y"
            mesh.attrib["face_node_connectivity"] = "Mesh2_face_nodes"
            mesh.attrib["face_dimension"]         = "nMesh2_face"

            # coordinates
            nx = defVar(ds, "Mesh2_node_x", Float64, ("nMesh2_node",))
            ny = defVar(ds, "Mesh2_node_y", Float64, ("nMesh2_node",))
            nx.attrib["standard_name"] = "longitude"
            nx.attrib["units"]         = "degrees_east"
            ny.attrib["standard_name"] = "latitude"
            ny.attrib["units"]         = "mdegrees_north"

            # connectivity
            FILL = Int32(2_147_483_647)
            f2n = defVar(ds, "Mesh2_face_nodes", Int32, ("nMesh2_face", "nMaxMesh2_face_nodes"); fillvalue=FILL)
            f2n.attrib["cf_role"]     = "face_node_connectivity"
            f2n.attrib["start_index"] = 1

            # data variables
            for ivar = 1:noutvar
                data_var = defVar(ds, "q$(ivar)", Float64, ("nMesh2_node",))
                data_var.attrib["long_name"] = "q$(ivar) field"
                data_var.attrib["location"]  = "node"
                data_var.attrib["mesh"]  = "mesh"
                data_var.attrib["units"]     = "N/A"

                data_var[:] = global_qout[:, ivar]
            end

            # global attributes
            ds.attrib["title"]       = "Unstructured data (MPI gathered)"
            ds.attrib["Conventions"] = "CF-1.11 UGRID-1.0"

            # write coordinate and connectivity data
            nx[:]  = global_xx
            ny[:]  = global_yy
            f2n[:] = global_subelem
            
        end
        
        println("Rank 0: Wrote global NetCDF file with $global_npoin nodes, $global_nsubelem elements")
    end
    
    MPI.Barrier(comm)
    
end

function write_NetCDF(SD::NSD_3D, mesh::St_mesh, q::Array, qaux::Array, mp, 
                      connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                      t, title::String, OUTPUT_DIR::String, inputs::Dict, varnames, outvarnames;
                      iout=1, nvar=1, qexact=zeros(1,nvar), case="")

    if (isa(varnames, Tuple)    || isa(varnames, String) )   varnames    = collect(varnames) end
    if (isa(outvarnames, Tuple) || isa(outvarnames, String)) outvarnames = collect(outvarnames) end
    
    nvar    = size(varnames, 1)
    noutvar = max(nvar, size(outvarnames,1))
    npoin   = mesh.npoin
    
    xx = zeros(size(mesh.x,1))
    yy = zeros(size(mesh.x,1))
    zz = zeros(size(mesh.x,1))
    xx .= mesh.x 
    yy .= mesh.y
    zz .= mesh.z
    nsubelem = mesh.nelem*(mesh.ngl-1)^3
    subelem  = Array{Int64}(undef, nsubelem, 8)
    
    isel = 1
    for iel = 1:mesh.nelem
        for i = 1:mesh.ngl-1
            for j = 1:mesh.ngl-1
                for k = 1:mesh.ngl-1
                    ip1 = mesh.connijk[iel,i,j,k]
                    ip2 = mesh.connijk[iel,i+1,j,k]
                    ip3 = mesh.connijk[iel,i+1,j+1,k]
                    ip4 = mesh.connijk[iel,i,j+1,k]
                    
                    ip5 = mesh.connijk[iel,i,j,k+1]
                    ip6 = mesh.connijk[iel,i+1,j,k+1]
                    ip7 = mesh.connijk[iel,i+1,j+1,k+1]
                    ip8 = mesh.connijk[iel,i,j+1,k+1]

                    subelem[isel, 1] = ip1
                    subelem[isel, 2] = ip2
                    subelem[isel, 3] = ip3
                    subelem[isel, 4] = ip4
                    subelem[isel, 5] = ip5
                    subelem[isel, 6] = ip6
                    subelem[isel, 7] = ip7
                    subelem[isel, 8] = ip8
                    
                    
                    isel = isel + 1
                end
            end
        end
    end
    
    #
    # Fetch user-defined diagnostic vars or take them from the solution vars:
    #
    qout = zeros(Float64, npoin, noutvar)
    @show noutvar
    u2uaux!(qaux, q, nvar, npoin)
    call_user_uout(qout, qaux, qexact, mp, inputs[:SOL_VARS_TYPE], npoin, nvar, noutvar)

    fout_name = string(OUTPUT_DIR, "/iter_", iout, ".nc")
    NCDataset(fout_name, "c") do ds
        
        # dimensions
        defDim(ds, "nMesh3d_node", npoin)
        defDim(ds, "nMesh3d_volume", nsubelem)
        defDim(ds, "nMaxMesh3d_volume_nodes", 8)

        # --- the mesh topology "dummy" variable with your exact attributes ---
        mesh = defVar(ds, "mesh", Int32, ())  # scalar dummy
        mesh.attrib["cf_role"]                  = "mesh_topology"
        mesh.attrib["long_name"]                = "Topology of a 3-d unstructured mesh"
        mesh.attrib["topology_dimension"]       = 3
        mesh.attrib["node_coordinates"]         = "Mesh3d_node_x Mesh3d_node_y Mesh3d_node_z"
        mesh.attrib["volume_node_connectivity"] = "Mesh3d_volume_nodes"
        mesh.attrib["volume_shape_type"]        = "Mesh3d_vol_types"
        mesh.attrib["volume_dimension"]         = "nMesh3d_volume"

        vol_types = defVar(ds, "Mesh3d_vol_types", Int32, ("nMesh3d_volume",))
        vol_types.attrib["cf_role"]       = "volume_shape_type"
        vol_types.attrib["long_name"]     = "Specifies the shape of the individual volumes."
        vol_types.attrib["flag_range"]    = [0, 2]                # integer array
        vol_types.attrib["flag_values"]   = [0, 1, 2]             # integer array
        vol_types.attrib["flag_meanings"] = "tetrahedron wedge hexahedron"


        # node coordinates
        nx = defVar(ds, "Mesh3d_node_x", Float64, ("nMesh3d_node",))
        ny = defVar(ds, "Mesh3d_node_y", Float64, ("nMesh3d_node",))
        nz = defVar(ds, "Mesh3d_node_z", Float64, ("nMesh3d_node",))
        nx.attrib["standard_name"] = "projection_x_coordinate"
        nx.attrib["units"]         = "m"
        ny.attrib["standard_name"] = "projection_y_coordinate"
        ny.attrib["units"]         = "m"
        nz.attrib["standard_name"] = "projection_z_coordinate"
        nz.attrib["units"]         = "m"

        # volum->node connectivity 
        FILL = Int32(2_147_483_647)
        v2n = defVar(ds, "Mesh3d_volume_nodes", Int32, ("nMesh3d_volume", "nMaxMesh3d_volume_nodes");fillvalue=FILL)
        v2n.attrib["cf_role"]     = "volume_node_connectivity"
        v2n.attrib["start_index"] = 1   # choose 1-based; UGRID default is 0-based if unspecified

        for ivar = 1:noutvar

            # Define data variable
            data_var = defVar(ds, "q$(ivar)", Float64, ("nMesh3d_node",))
            
            # Add attributes for data
            data_var.attrib["long_name"] = "q$(ivar) field"
            data_var.attrib["location"]  = "node"
            data_var.attrib["units"]     = "N/A"
            data_var[:]                  = qout[:,ivar]
        end
            
        # Add global attributes
        ds.attrib["title"] = "Unstructured data"
        ds.attrib["Conventions"] = "CF-1.11 UGRID-1.0"

        volume_type = zeros(Int32, nsubelem)
        fill!(volume_type, 2)
        vol_types[:] = volume_type
        nx[:]        = xx
        ny[:]        = yy
        nz[:]        = zz
        v2n[:]       = subelem
            
    end

end

#---------------------------------------------------------------------------------------------------------------------
# Write structured VTK on Alya grid
#---------------------------------------------------------------------------------------------------------------------
using Base64

# ═══════════════════════════════════════════════════════════════════════════ #
#  Entry point                                                                #
# ═══════════════════════════════════════════════════════════════════════════ #
function write_pvtu_structured(
    filename_base::String,
    alya_coords::Matrix{Float64},
    u_interp::Matrix{Float64},
    field_names::Vector{String};
    nparts::Int   = 1,
    time::Float64 = 0.0,
)
    n_points, neqs = size(u_interp)
    @assert size(alya_coords, 1) == n_points "coords/u_interp row mismatch: $(size(alya_coords,1)) vs $n_points"
    @assert 1 <= length(field_names) <= neqs  "need 1..neqs field names"

    nx, ny = _detect_grid_dims(alya_coords)
    @info "Detected structured grid: nx=$nx  ny=$ny  (total=$(nx*ny), actual=$n_points)"
    if nx * ny != n_points
        error("Grid detection failed: nx=$nx × ny=$ny = $(nx*ny) ≠ $n_points points.\n" *
              "x range: [$(minimum(alya_coords[:,1])), $(maximum(alya_coords[:,1]))]\n" *
              "y range: [$(minimum(alya_coords[:,2])), $(maximum(alya_coords[:,2]))]")
    end

    perm = _structured_permutation(alya_coords, nx, ny)

    dir  = dirname(filename_base)
    isempty(dir) || mkpath(dir)
    outdir = isempty(dir) ? "." : dir
    base   = basename(filename_base)

    # Write each VTU piece
    vtu_names = String[]
    for ip in 1:nparts
        jrange   = _part_jrange(ny, nparts, ip)
        vtu_name = "$(base)_$(lpad(ip, 4, '0')).vtu"
        push!(vtu_names, vtu_name)
        vtu_path = joinpath(outdir, vtu_name)
        xml = _build_vtu_xml(alya_coords, u_interp, field_names, perm, nx, ny, jrange)
        write(vtu_path, xml)          # atomic: write string, not stream
        @info "Wrote VTU piece $ip/$nparts: $vtu_path"
    end

    # Write PVTU master
    pvtu_path = filename_base * ".pvtu"
    xml = _build_pvtu_xml(vtu_names, field_names, time)
    write(pvtu_path, xml)             # atomic: write string, not stream
    @info "Wrote PVTU: $pvtu_path"
    return pvtu_path
end


# ═══════════════════════════════════════════════════════════════════════════ #
#  XML builders — return complete Strings, never use incremental IO          #
# ═══════════════════════════════════════════════════════════════════════════ #
function _build_pvtu_xml(
    vtu_names::Vector{String},
    field_names::Vector{String},
    time::Float64,
)
    io = IOBuffer()
    println(io, """<?xml version="1.0"?>""")
    println(io, """<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">""")
    println(io, """  <PUnstructuredGrid GhostLevel="0">""")
    println(io, """    <PFieldData>""")
    println(io, """      <PDataArray type="Float64" Name="TimeValue" NumberOfTuples="1"/>""")
    println(io, """    </PFieldData>""")
    println(io, """    <PPointData>""")
    for name in field_names
        println(io, """      <PDataArray type="Float64" Name="$name"/>""")
    end
    println(io, """    </PPointData>""")
    println(io, """    <PPoints>""")
    println(io, """      <PDataArray type="Float64" NumberOfComponents="3"/>""")
    println(io, """    </PPoints>""")
    # ← this is the critical block that was ending up after </VTKFile>
    for vtu in vtu_names
        println(io, """    <Piece Source="$(basename(vtu))"/>""")
    end
    println(io, """  </PUnstructuredGrid>""")
    println(io, """</VTKFile>""")
    return String(take!(io))
end


function _build_vtu_xml(
    alya_coords::Matrix{Float64},
    u_interp::Matrix{Float64},
    field_names::Vector{String},
    perm::Vector{Int},
    nx::Int, ny::Int,
    jrange::UnitRange{Int},
)
    local_ny    = length(jrange)
    n_pts       = nx * local_ny
    n_cells     = (nx - 1) * (local_ny - 1)

    # Local permutation: structured (i,j_local) → global alya_coords row
    local_perm = Vector{Int}(undef, n_pts)
    for (lj, j) in enumerate(jrange), i in 1:nx
        local_perm[(lj-1)*nx + i] = perm[(j-1)*nx + i]
    end

    # 3-component coordinates
    pts = Vector{Float64}(undef, 3 * n_pts)
    for lk in 1:n_pts
        gk = local_perm[lk]
        pts[3*(lk-1)+1] = alya_coords[gk, 1]
        pts[3*(lk-1)+2] = alya_coords[gk, 2]
        pts[3*(lk-1)+3] = 0.0
    end

    # Quad connectivity (0-based VTK indices)
    connectivity = Vector{Int32}(undef, 4 * n_cells)
    offsets      = Vector{Int32}(undef, n_cells)
    types        = fill(UInt8(9), n_cells)   # VTK_QUAD
    for (c, (jl, i)) in enumerate(Iterators.product(1:(local_ny-1), 1:(nx-1)))
        p00 = Int32((jl-1)*nx + (i-1))
        p10 = Int32((jl-1)*nx +  i   )
        p11 = Int32( jl   *nx +  i   )
        p01 = Int32( jl   *nx + (i-1))
        connectivity[4*(c-1)+1] = p00
        connectivity[4*(c-1)+2] = p10
        connectivity[4*(c-1)+3] = p11
        connectivity[4*(c-1)+4] = p01
        offsets[c] = Int32(4*c)
    end

    io = IOBuffer()
    println(io, """<?xml version="1.0"?>""")
    println(io, """<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt64">""")
    println(io, """  <UnstructuredGrid>""")
    println(io, """    <Piece NumberOfPoints="$n_pts" NumberOfCells="$n_cells">""")

    println(io, """      <Points>""")
    println(io, """        <DataArray type="Float64" NumberOfComponents="3" format="binary">""")
    println(io, "          ", _b64(pts))
    println(io, """        </DataArray>""")
    println(io, """      </Points>""")

    println(io, """      <PointData>""")
    for (q, name) in enumerate(field_names)
        arr = Float64[u_interp[local_perm[lk], q] for lk in 1:n_pts]
        println(io, """        <DataArray type="Float64" Name="$name" format="binary">""")
        println(io, "          ", _b64(arr))
        println(io, """        </DataArray>""")
    end
    println(io, """      </PointData>""")

    println(io, """      <Cells>""")
    println(io, """        <DataArray type="Int32" Name="connectivity" format="binary">""")
    println(io, "          ", _b64(connectivity))
    println(io, """        </DataArray>""")
    println(io, """        <DataArray type="Int32" Name="offsets" format="binary">""")
    println(io, "          ", _b64(offsets))
    println(io, """        </DataArray>""")
    println(io, """        <DataArray type="UInt8" Name="types" format="binary">""")
    println(io, "          ", _b64(types))
    println(io, """        </DataArray>""")
    println(io, """      </Cells>""")

    println(io, """    </Piece>""")
    println(io, """  </UnstructuredGrid>""")
    println(io, """</VTKFile>""")
    return String(take!(io))
end


# ═══════════════════════════════════════════════════════════════════════════ #
#  Binary base64 block: UInt64 byte-count header + raw data                  #
# ═══════════════════════════════════════════════════════════════════════════ #
function _b64(data::AbstractVector{T}) where {T}
    buf = IOBuffer()
    write(buf, UInt64(length(data) * sizeof(T)))
    for x in data
        write(buf, x)
    end
    return base64encode(take!(buf))
end


# ═══════════════════════════════════════════════════════════════════════════ #
#  Grid helpers                                                               #
# ═══════════════════════════════════════════════════════════════════════════ #
function _detect_grid_dims(coords::Matrix{Float64})
    function find_n_unique(vals)
        sv = sort(vals)
        # Find the smallest nonzero gap between adjacent values
        min_gap = Inf
        for k in 2:length(sv)
            d = sv[k] - sv[k-1]
            d > 0 && (min_gap = min(min_gap, d))
        end
        # Tolerance = half the minimum spacing — robust regardless of domain size
        atol = min_gap / 2
        n = 1
        for k in 2:length(sv)
            sv[k] - sv[k-1] > atol && (n += 1)
        end
        return n, atol
    end

    nx, atol_x = find_n_unique(coords[:, 1])
    ny, atol_y = find_n_unique(coords[:, 2])
    @info "Grid detection: nx=$nx (atol_x=$atol_x)  ny=$ny (atol_y=$atol_y)"
    return nx, ny
end

function _cluster_values(vals::Vector{Float64}, n::Int)
    sv = sort(vals)
    min_gap = Inf
    for k in 2:length(sv)
        d = sv[k] - sv[k-1]
        d > 0 && (min_gap = min(min_gap, d))
    end
    atol = min_gap / 2

    reps = Float64[sv[1]]
    for v in sv[2:end]
        v - reps[end] > atol ? push!(reps, v) : (reps[end] = (reps[end] + v) / 2)
    end
    length(reps) == n || error("clustering gave $(length(reps)) groups, expected $n")
    return reps
end

function _nearest_grid_idx(sv::Vector{Float64}, v::Float64)
    idx = searchsortedfirst(sv, v)
    idx > length(sv) && (idx = length(sv))
    idx > 1 && abs(sv[idx-1] - v) < abs(sv[idx] - v) && (idx -= 1)
    return idx
end

function _structured_permutation(coords::Matrix{Float64}, nx::Int, ny::Int)
    xs = coords[:, 1]
    ys = coords[:, 2]
    x_uniq = _cluster_values(xs, nx)
    y_uniq = _cluster_values(ys, ny)

    perm = zeros(Int, nx * ny)
    for k in 1:size(coords, 1)
        i   = _nearest_grid_idx(x_uniq, xs[k])
        j   = _nearest_grid_idx(y_uniq, ys[k])
        lin = (j-1)*nx + i
        perm[lin] != 0 && @warn "Slot ($i,$j) already filled — duplicate coordinate at row $k"
        perm[lin] = k
    end
    any(==(0), perm) && error("$(count(==(0),perm)) grid slots unfilled after permutation")
    return perm
end

function _part_jrange(ny::Int, nparts::Int, ip::Int)
    base  = ny ÷ nparts
    extra = ny % nparts
    j0 = (ip-1)*base + min(ip-1, extra) + 1
    j1 = j0 + base - 1 + (ip <= extra ? 1 : 0)
    return j0:j1
end
