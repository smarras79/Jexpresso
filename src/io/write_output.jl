using WriteVTK

include("./plotting/jeplots.jl")

#
# VTK 2D/3D
#
function write_output(SD, sol::ODESolution, mesh::St_mesh, mp, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::VTK; nvar=1, qexact=zeros(1,nvar), case="")
 
    for iout = 1:size(sol.t[:],1)
        if (inputs[:backend] == CPU())
            title = @sprintf "final solution at t=%6.4f" sol.t[iout]
            write_vtk(SD, mesh, sol.u[iout][:], mp, sol.t[iout], title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, qexact=qexact, case=case)
        else
            u = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin*(nvar+1))
            KernelAbstractions.copyto!(CPU(),u,sol.u[iout][:])
            u_exact = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin,nvar+1)
            KernelAbstractions.copyto!(CPU(),u_exact,qexact)
            convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
            title = @sprintf "final solution at t=%6.4f" sol.t[iout]
            write_vtk(SD, mesh, u, mp, sol.t[iout], title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, qexact=u_exact, case=case)
        end
    end
    println(string(" # Writing output to VTK file:", OUTPUT_DIR, "*.vtu ... DONE") )
end

function write_output(SD, sol::SciMLBase.LinearSolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::VTK; nvar=1, qexact=zeros(1,nvar), case="")
    #
    # 2D VTK of x from Ax=b
    #    
    if (inputs[:backend] == CPU())
        title = @sprintf "Solution"
        write_vtk(SD, mesh, sol.u, "1", 0, title, OUTPUT_DIR, inputs, varnames; iout=1, nvar=nvar, qexact=qexact, case="")
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

function write_output(SD, u_sol, t, iout, mesh::St_mesh, mp, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::VTK; nvar=1, qexact=zeros(1,nvar), case="")
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    title = @sprintf "final solution at t=%6.4f" iout
    if (inputs[:backend] == CPU())
        write_vtk(SD, mesh, u_sol, mp, t, title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, qexact=qexact, case=case)        
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

function write_vtk(SD::NSD_2D, mesh::St_mesh, q::Array, mp, t, title::String, OUTPUT_DIR::String, inputs::Dict, varnames; iout=1, nvar=1, qexact=zeros(1,nvar), case="")

    outvars = varnames
    nvars   = length(outvars)
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

    qaux = zeros(Float64, npoin, nvar+1)
    qout = zeros(Float64, npoin, nvar+1)
    u2uaux!(qaux, q, nvar, npoin)
    for ip=1:npoin
        user_uout!(@view(qout[ip,1:nvar]), @view(qaux[ip,1:nvar]), qexact[ip,1:nvar], inputs[:SOL_VARS_TYPE])
    end    
    uaux2u!(q, qout, nvar, npoin)
    
    #
    # Write solution:
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
        
        for ivar = 1:nvar
            idx = (ivar - 1)*npoin
            vtkf[string(varnames[ivar]), VTKPointData()] = @view(q[idx+1:ivar*npoin])
        end
        vtkf
    end
    
    outfiles = map(vtk_save, vtkfile)
    
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
    
    vtkfile = vtk_grid(fout_name, mesh.x[1:mesh.npoin], mesh.y[1:mesh.npoin], mesh.y[1:mesh.npoin]*TFloat(0.0), cells)

    for ivar = 1:length(outvarsref)
        vtkfile[string(outvarsref[ivar]), VTKPointData()] =  @view(q[1:mesh.npoin,ivar])
    end
    outfiles = vtk_save(vtkfile)


    vtkfile = map(mesh.parts) do part
        vtkf = pvtk_grid(file_name, mesh.x[1:mesh.npoin], mesh.y[1:mesh.npoin], mesh.z[1:mesh.npoin], cells, compress=false;
                        part=part, nparts=nparts, ismain=(part==1))
        vtkf["part", VTKCellData()] = ones(isel -1) * part

        for ivar = 1:length(outvarsref)
            vtkf[string(outvarsref[ivar]), VTKPointData()] =  @view(q[1:mesh.npoin,ivar])
        end
        vtkf
    end
    outfiles = map(vtk_save, vtkfile)
end



function write_vtk_ref(SD::NSD_3D, mesh::St_mesh, q::Array, file_name::String, OUTPUT_DIR::String; iout=1, nvar=1, qexact=zeros(1,nvar), case="", outvarsref=tuple(("" for _ in 1:nvar)))

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
    
    vtkfile = vtk_grid(fout_name, mesh.x[1:mesh.npoin], mesh.y[1:mesh.npoin], mesh.z[1:mesh.npoin], cells)

    for ivar = 1:length(outvarsref)
        vtkfile[string(outvarsref[ivar]), VTKPointData()] =  @view(q[1:mesh.npoin,ivar])
    end
    outfiles = vtk_save(vtkfile)


    vtkfile = map(mesh.parts) do part
        vtkf = pvtk_grid(file_name, mesh.x[1:mesh.npoin], mesh.y[1:mesh.npoin], mesh.z[1:mesh.npoin], cells, compress=false;
                        part=part, nparts=nparts, ismain=(part==1))
        vtkf["part", VTKCellData()] = ones(isel -1) * part

        for ivar = 1:length(outvarsref)
            vtkf[string(outvarsref[ivar]), VTKPointData()] =  @view(q[1:mesh.npoin,ivar])
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
