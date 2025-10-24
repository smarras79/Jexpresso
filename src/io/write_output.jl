using WriteVTK

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
        fig = Figure(size = (1200,800),fontsize=22)
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
        fig = Figure(size = (1200,800),fontsize=22)
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

    comm = MPI.COMM_WORLD
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
                      connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                      OUTPUT_DIR::String, inputs::Dict,
                      varnames, outvarnames,
                      outformat::VTK;
                      nvar=1, qexact=zeros(1,nvar), case="")
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    title = @sprintf "final solution at t=%6.4f" iout
    if (inputs[:backend] == CPU())

        write_vtk(SD, mesh, sol, uaux, mp, 
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

    comm = MPI.COMM_WORLD
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
                   connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                   t, title::String, OUTPUT_DIR::String, inputs::Dict, varnames, outvarnames;
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
    
    comm = MPI.COMM_WORLD
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
    comm = MPI.COMM_WORLD
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


    comm = MPI.COMM_WORLD
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