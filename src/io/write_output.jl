using WriteVTK

include("./plotting/jeplots.jl")

#------------------------------------------------------------------
# Callback for missing user_uout!()
#------------------------------------------------------------------
function call_user_uout(uout, u, qe, mp, ET, npoin, nvar, noutvar)
    
    if function_exists(@__MODULE__, :user_uout!)
        for ip=1:npoin
            user_uout!(ip, ET, @view(uout[ip,1:noutvar]), @view(u[ip,1:nvar]), @view(qe[ip,1:nvar]), mp)
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
    noutvar  = max(nvar, size(outvarnames,1))
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
    noutvar = max(nvar, size(outvarnames,1))
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
function write_output(SD, u::AbstractArray, t, iout, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::HDF5; nvar=1, qexact=zeros(1,nvar), case="")
    
    # println(string(" # Writing restart HDF5 file:", OUTPUT_DIR, "*.h5 ...  ") )
    iout = size(t,1)
    title = @sprintf "Final solution at t=%6.4f" t
    if (inputs[:backend] == CPU())
    
        write_hdf5(SD, mesh, u, qexact, title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, case=case)
    else
        u_gpu = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin*nvar)
        KernelAbstractions.copyto!(CPU(),u_gpu, u)
        u_exact = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin,nvar+1)
        KernelAbstractions.copyto!(CPU(),u_exact,qexact)
        convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
        write_hdf5(SD, mesh, u_gpu, u_exact, title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, case=case)
        convert_mesh_arrays!(SD, mesh, inputs[:backend], inputs)
    end

    
    println(string(" # Writing restart HDF5 file:", OUTPUT_DIR, "*.h5 ... DONE") )
    
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
function read_output(SD::NSD_2D, INPUT_DIR::String, inputs::Dict, npoin, outformat::HDF5; nvar=1)
    
    #println(string(" # Reading restart HDF5 file:", INPUT_DIR, "*.h5 ...  ") )
    q, qe = read_hdf5(SD, INPUT_DIR, inputs, npoin, nvar)
    println(string(" # Reading restart HDF5 file:", INPUT_DIR, "*.h5 ... DONE") )

    return q, qe
end


function write_hdf5(SD, mesh::St_mesh, q::AbstractArray, qe::AbstractArray, title::String, OUTPUT_DIR::String, inputs::Dict, varnames; iout=1, nvar=1, case="")
    
    #Write one HDF5 file per variable
    for ivar = 1:nvar
        fout_name = string(OUTPUT_DIR, "/var_", ivar, ".h5")
        idx = (ivar - 1)*mesh.npoin
        
        h5open(fout_name, "w") do fid        
            write(fid, "q",  q[idx+1:ivar*mesh.npoin]);
            write(fid, "qe", qe[1:mesh.npoin, ivar]);
        end

    end
end

function read_hdf5(SD, INPUT_DIR::String, inputs::Dict, npoin, nvar)
    
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
