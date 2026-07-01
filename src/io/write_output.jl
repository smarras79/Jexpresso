using WriteVTK

include("./plotting/jeplots.jl")

#------------------------------------------------------------------
# Write PVD file for ParaView time series
#------------------------------------------------------------------
function init_pvd_file(path)
    open(path, "w") do io
        println(io, "<?xml version=\"1.0\"?>")
        println(io, "<VTKFile type=\"Collection\" version=\"0.1\">")
        println(io, "  <Collection>")
        println(io, "  </Collection>")
        println(io, "</VTKFile>")
    end
end

function append_pvd_entry(path, time, filename)
    lines = readlines(path)
    insert_pos = length(lines) - 1  # insert before last 2 lines
    insert!(lines, insert_pos, "    <DataSet timestep=\"$time\" file=\"$filename\"/>")
    open(path, "w") do io
        for line in lines
            println(io, line)
        end
    end
end

#------------------------------------------------------------------
# Callback for missing user_uout!()
#------------------------------------------------------------------
function call_user_uout(uout, u, qe, mp, ET, npoin, nvar, noutvar; μ_dsgs_pnode=nothing)

    if function_exists(@__MODULE__, :user_uout!)
        for ip=1:npoin
            user_uout!(ip, ET, @view(uout[ip,1:noutvar]), @view(u[ip,:]), @view(qe[ip,:]);
                       mp=mp, μ_dsgs_pnode=μ_dsgs_pnode)
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

function write_output(SD::NSD_1D, q::Array, t, iout, mesh::St_mesh, OUTPUT_DIR::String, inputs, varnames, outformat::PNG; nvar=1, qexact=zeros(1,nvar), case="")
    #OK
    nvar = length(varnames)
    qout = zeros(mesh.npoin)
    
    plot_results(SD, mesh, q[:], "initial", OUTPUT_DIR, varnames, inputs; iout=1, nvar=nvar, PT=nothing)
end

function write_output(SD::NSD_1D, sol, uaux, t, iout,  mesh::St_mesh, mp,
                      connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                      OUTPUT_DIR::String, inputs,
                      varnames, outvarnames,
                      outformat::PNG;
                      nvar=1, qexact=zeros(1,nvar), case="",
                      μ_dsgs_pnode=nothing)
        
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
        title = @sprintf "t = %.4f" t
        # DSGS runs render the viscosity staircase as one more panel of
        # the same output time (the per-node broadcast is in μ_dsgs_pnode)
        μ_nodes = (μ_dsgs_pnode !== nothing && inputs[:backend] == CPU()) ? μ_dsgs_pnode : nothing
            if (inputs[:backend] == CPU())
                plot_results(SD, mesh, sol, title, OUTPUT_DIR, varnames, inputs; iout=iout, nvar=nvar, PT=nothing, μ_nodes=μ_nodes)
            else
                uout = KernelAbstractions.allocate(CPU(), TFloat, Int64(mesh.npoin*nvar))
                KernelAbstractions.copyto!(CPU(), uout, sol)
                convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
                plot_results(SD, mesh, uout, title, OUTPUT_DIR, varnames, inputs; iout=iout, nvar=nvar, PT=nothing, μ_nodes=μ_nodes)
            end
        #end
    end
    MPI.Comm_rank(get_mpi_comm()) == 0 && println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE ") )
end

function write_output(SD::NSD_2D, sol, uaux, t, iout,  mesh::St_mesh, mp,
                      connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                      OUTPUT_DIR::String, inputs,
                      varnames, outvarnames,
                      outformat::PNG;
                      nvar=1, qexact=zeros(1,nvar), case="",
                      μ_dsgs_pnode=nothing)

    #
    # 2D PNG of q(t): one colored map per variable and output time.
    # inputs[:lplot_surf3d] selects the Spline2D-interpolated surface view
    # (plot_surf3d), otherwise the nodal point map (plot_triangulation) is
    # used -- the latter is the safer choice for solutions with kinks such
    # as the shallow water wet/dry front, where a global spline overshoots.
    #
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)

    if (inputs[:backend] == CPU())
        q = sol
    else
        q = KernelAbstractions.allocate(CPU(), TFloat, Int64(mesh.npoin*nvar))
        KernelAbstractions.copyto!(CPU(), q, sol)
        convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
    end

    title = @sprintf "t = %.4f s" t
    if (inputs[:lplot_surf3d])
        plot_surf3d(SD, mesh, q, title, OUTPUT_DIR;
                    iout=iout, nvar=nvar,
                    smoothing_factor=inputs[:smoothing_factor], varnames=varnames)
    else
        plot_triangulation(SD, mesh, q, title, OUTPUT_DIR, inputs;
                           iout=iout, nvar=nvar, varnames=varnames)
    end

    println_rank(string(" # writing ", OUTPUT_DIR, "/<var>-it", iout, ".png at t=", t, " s... DONE"); msg_rank = rank)
end


function write_output(SD, sol::SciMLBase.LinearSolution, uaux, mesh::St_mesh,
                      OUTPUT_DIR::String, inputs,
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
    
    MPI.Comm_rank(get_mpi_comm()) == 0 && println(string(" # Writing output to VTK file:", OUTPUT_DIR, "*.pvtu ... DONE") )

end


function write_output(SD, sol, uaux, t, iout,  mesh::St_mesh, mp,
                      connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                      OUTPUT_DIR::String, inputs,
                      varnames, outvarnames,
                      outformat::VTK;
                      nvar=1, qexact=zeros(1,nvar), case="",
                      μ_dsgs_pnode=nothing, metrics=nothing,
                      extra_fields=Pair{String,Vector{Float64}}[])

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    title = @sprintf "final solution at t=%6.4f" iout
    if (inputs[:backend] == CPU())

        write_vtk(SD, mesh, sol, uaux, mp,
                  connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                  t, title, OUTPUT_DIR, inputs,
                  varnames, outvarnames;
                  iout=iout, nvar=nvar, qexact=qexact, case=case,
                  μ_dsgs_pnode=μ_dsgs_pnode, metrics=metrics,
                  extra_fields=extra_fields)
        
    else
        #VERIFY THIS on GPU
        u = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin*(nvar+1))
        KernelAbstractions.copyto!(CPU(),u, u_sol)
        u_exact = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin,nvar+1)
        KernelAbstractions.copyto!(CPU(),u_exact,qexact)
        convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
        write_vtk(SD, mesh, u, mp, t, title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, qexact=u_exact, case=case)
    end

    println_rank(string(" # writing ", OUTPUT_DIR, "/iter_", iout, ".pvtu at t=", t, " s... DONE"); msg_rank = rank )

end

function write_output(SD, sol, uaux, t, iout,  mesh::St_mesh, mp,
                    connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                    OUTPUT_DIR::String, inputs,
                    varnames, outvarnames,
                    outformat::NETCDF;
                    nvar=1, qexact=zeros(1,nvar), case="",
                    μ_dsgs_pnode=nothing)

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

    println_rank(string(" # writing ", OUTPUT_DIR, "/iter_", iout, ".pvtu at t=", t, " s... DONE"); msg_rank = rank )

end

#------------
# VTK writer
#------------
function write_vtk(SD::NSD_2D, mesh::St_mesh, q::Array, qaux::Array, mp,
                   connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                   t, title::String, OUTPUT_DIR::String, inputs, varnames, outvarnames;
                   iout=1, nvar=1, qexact=zeros(1,nvar), case="",
                   μ_dsgs_pnode=nothing, metrics=nothing,
                   extra_fields=Pair{String,Vector{Float64}}[])

    if (isa(varnames, Tuple)    || isa(varnames, String) )   varnames    = collect(varnames) end
    if (isa(outvarnames, Tuple) || isa(outvarnames, String)) outvarnames = collect(outvarnames) end
    
    nvar     = size(varnames, 1)
    noutvar  = size(outvarnames,1) #max(nvar, size(outvarnames,1))
    new_size = size(mesh.x,1)

    npoin          = mesh.npoin
    nelem          = mesh.nelem
    nelem_semi_inf = mesh.nelem_semi_inf
    ngl            = mesh.ngl
    ngr            = mesh.ngr
    
    if (nelem_semi_inf > 0)
        subelem = Array{Int64}(undef, nelem*(ngl-1)^2+nelem_semi_inf*(ngl-1)*(ngr-1), 4)
        cells = [MeshCell(VTKCellTypes.VTK_QUAD, [1, 2, 4, 3]) for _ in 1:nelem*(ngl-1)^2+nelem_semi_inf*(ngl-1)*(ngr-1)]
    else
        subelem = Array{Int64}(undef, nelem*(ngl-1)^2, 4)
        cells = [MeshCell(VTKCellTypes.VTK_QUAD, [1, 2, 4, 3]) for _ in 1:mesh.nelem*(ngl-1)^2]
    end
    
    isel = 1
    for iel = 1:nelem
        for i = 1:ngl-1
            for j = 1:ngl-1
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

    for iel = 1:nelem_semi_inf
        for i = 1:ngl-1
            for j = 1:ngr-1
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
    call_user_uout(qout, qaux, qexact, mp, inputs[:SOL_VARS_TYPE], npoin, nvar, noutvar;
                   μ_dsgs_pnode=μ_dsgs_pnode)


    #
    # Optional: non-constant diffusivity field a(x,y), written to the VTU so the
    # variable coefficient can be visualised alongside the solution. Active only
    # when inputs[:diffusivity] is a callable (x,y)->a (e.g. case1_nonconstant);
    # other cases are unaffected.
    #
    write_diffusivity = haskey(inputs, :diffusivity) && inputs[:diffusivity] !== nothing
    diffvals = write_diffusivity ?
        Float64[inputs[:diffusivity](mesh.x[ip], mesh.y[ip]) for ip = 1:npoin] :
        Float64[]

    #
    # Optional: GEOMETRY-induced reference diffusivity â (EL notes), written to
    # the VTU so the non-constant diffusion dictated by the element shapes can be
    # visualised. â(ξ) = (|K|/|ref|)·a·J_K⁻¹J_K⁻ᵀ with physical a = 1, so
    # det(â) = 1 (area-normalised) and the variation is purely anisotropic.
    #
    # Active when inputs[:lEL_nonconstant] is true AND metrics are supplied.
    # Three output formats, selected by inputs[:ahat_output] (default :cell):
    #   :cell   → per-cell invariants  ahat_eff = tr(â)/2 ≥ 1,  ahat_aniso = λmax/λmin
    #             (element-averaged; one value per element, piecewise-constant).
    #   :nodal  → the SAME invariants as a SMOOTHED nodal field (â averaged over the
    #             elements sharing each node, then ahat_eff / ahat_aniso from it).
    #   :tensor → the full smoothed nodal tensor components ahat_11, ahat_12, ahat_22.
    #
    write_ahat = get(inputs, :lEL_nonconstant, false) && metrics !== nothing
    ahat_mode  = get(inputs, :ahat_output, :cell)

    ahat_eff_cell   = Float64[];  ahat_aniso_cell = Float64[]   # :cell
    ahat_eff_node   = Float64[];  ahat_aniso_node = Float64[]   # :nodal
    ahat11_node     = Float64[];  ahat12_node     = Float64[];  ahat22_node = Float64[]  # :tensor

    if write_ahat
        if ahat_mode == :cell
            # ── per-cell invariants (element-averaged) ────────────────────────
            ncells = isel - 1
            ahat_eff_cell   = zeros(Float64, ncells)
            ahat_aniso_cell = zeros(Float64, ncells)
            ic = 1
            @inbounds for iel = 1:nelem
                seff = 0.0; saniso = 0.0
                for i = 1:ngl, j = 1:ngl
                    ix = metrics.dξdx[iel,i,j]; iy = metrics.dξdy[iel,i,j]
                    nx = metrics.dηdx[iel,i,j]; ny = metrics.dηdy[iel,i,j]
                    Je = metrics.Je[iel,i,j]
                    a11 = Je*(ix*ix + iy*iy); a12 = Je*(ix*nx + iy*ny); a22 = Je*(nx*nx + ny*ny)
                    tr   = a11 + a22
                    detâ = a11*a22 - a12*a12
                    disc = sqrt(max(tr*tr - 4*detâ, 0.0))
                    λmax = 0.5*(tr + disc); λmin = 0.5*(tr - disc)
                    seff   += 0.5*tr
                    saniso += (λmin > 0 ? λmax/λmin : 0.0)
                end
                eff = seff/(ngl*ngl); aniso = saniso/(ngl*ngl)
                for _ = 1:(ngl-1)*(ngl-1)
                    ahat_eff_cell[ic] = eff; ahat_aniso_cell[ic] = aniso; ic += 1
                end
            end
        else
            # ── smoothed NODAL â: average the tensor over elements at each node ─
            s11 = zeros(Float64, npoin); s12 = zeros(Float64, npoin)
            s22 = zeros(Float64, npoin); cnt = zeros(Int, npoin)
            @inbounds for iel = 1:nelem, i = 1:ngl, j = 1:ngl
                ip = mesh.connijk[iel,i,j]
                ix = metrics.dξdx[iel,i,j]; iy = metrics.dξdy[iel,i,j]
                nx = metrics.dηdx[iel,i,j]; ny = metrics.dηdy[iel,i,j]
                Je = metrics.Je[iel,i,j]
                s11[ip] += Je*(ix*ix + iy*iy)
                s12[ip] += Je*(ix*nx + iy*ny)
                s22[ip] += Je*(nx*nx + ny*ny)
                cnt[ip] += 1
            end
            @inbounds for ip = 1:npoin
                c = max(cnt[ip], 1); s11[ip] /= c; s12[ip] /= c; s22[ip] /= c
            end
            if ahat_mode == :tensor
                ahat11_node = s11; ahat12_node = s12; ahat22_node = s22
            else  # :nodal invariants from the smoothed tensor
                ahat_eff_node   = zeros(Float64, npoin)
                ahat_aniso_node = zeros(Float64, npoin)
                @inbounds for ip = 1:npoin
                    tr   = s11[ip] + s22[ip]
                    detâ = s11[ip]*s22[ip] - s12[ip]*s12[ip]
                    disc = sqrt(max(tr*tr - 4*detâ, 0.0))
                    λmax = 0.5*(tr + disc); λmin = 0.5*(tr - disc)
                    ahat_eff_node[ip]   = 0.5*tr
                    ahat_aniso_node[ip] = (λmin > 0 ? λmax/λmin : 0.0)
                end
            end
        end
    end

    #
    # Write solution to vtk:
    #
    fout_name = string(OUTPUT_DIR, "/iter_", iout)
    vtkfile = map(mesh.parts) do part
        vtkf = pvtk_grid(fout_name,
                         mesh.coords[1:mesh.npoin,1],
                         mesh.coords[1:mesh.npoin,2],
                         mesh.coords[1:mesh.npoin,2]*TFloat(0.0),
                         cells,
                         compress=false;
                         part=part, nparts=mesh.nparts, ismain=(part==1))
        vtkf["part", VTKCellData()] = ones(isel -1) * part

        for ivar = 1:noutvar
            idx = (ivar - 1)*npoin
            vtkf[string(outvarnames[ivar]), VTKPointData()] = @view(qout[1:npoin,ivar])
        end

        if write_diffusivity
            vtkf["diffusivity", VTKPointData()] = @view(diffvals[1:npoin])
        end

        if write_ahat
            if ahat_mode == :cell
                vtkf["ahat_eff",   VTKCellData()] = ahat_eff_cell
                vtkf["ahat_aniso", VTKCellData()] = ahat_aniso_cell
            elseif ahat_mode == :tensor
                vtkf["ahat_11", VTKPointData()] = @view(ahat11_node[1:npoin])
                vtkf["ahat_12", VTKPointData()] = @view(ahat12_node[1:npoin])
                vtkf["ahat_22", VTKPointData()] = @view(ahat22_node[1:npoin])
            else # :nodal
                vtkf["ahat_eff",   VTKPointData()] = @view(ahat_eff_node[1:npoin])
                vtkf["ahat_aniso", VTKPointData()] = @view(ahat_aniso_node[1:npoin])
            end
        end

        # Optional caller-supplied nodal fields (name => length-npoin vector).
        # Used by element learning to write the reference (direct SEM /
        # manufactured) solution and its difference from the inferred solution
        # alongside the solution itself, for side-by-side visualisation.
        for (fname, fdata) in extra_fields
            vtkf[fname, VTKPointData()] = @view(fdata[1:npoin])
        end

        vtkf
    end
    
    outfiles = map(vtk_save, vtkfile)
    
end

function write_vtk(SD::NSD_3D, mesh::St_mesh, q::Array, qaux::Array, mp,
                   connijk_original, poin_in_bdy_face_original,
                   x_original, y_original, z_original,
                   t, title::String, OUTPUT_DIR::String, inputs,
                   varnames, outvarnames;
                   iout=1, nvar=1, qexact=zeros(1,nvar), case="",
                   μ_dsgs_pnode=nothing,
                   extra_fields=Pair{String,Vector{Float64}}[])

    if (isa(varnames, Tuple)    || isa(varnames, String) )   varnames    = collect(varnames) end
    if (isa(outvarnames, Tuple) || isa(outvarnames, String)) outvarnames = collect(outvarnames) end
    
    nvar    = size(varnames, 1)
    noutvar = size(outvarnames,1) #max(nvar, size(outvarnames,1))
    npoin   = mesh.npoin
    
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
    call_user_uout(qout, qaux, qexact, mp, inputs[:SOL_VARS_TYPE], npoin, nvar, noutvar;
                   μ_dsgs_pnode=μ_dsgs_pnode)


    #
    # Write solution:
    #
    fout_name = string(OUTPUT_DIR, "/iter_", iout)
    vtkfile = map(mesh.parts) do part
        vtkf = pvtk_grid(fout_name,
                         mesh.coords[1:mesh.npoin,1],
                         mesh.coords[1:mesh.npoin,2],
                         mesh.coords[1:mesh.npoin,3],
                         cells,
                         compress=false;
                         part=part, nparts=mesh.nparts, ismain=(part==1))
        vtkf["part", VTKCellData()] = ones(isel -1) * part

        for ivar = 1:noutvar
            idx = (ivar - 1)*npoin
            vtkf[string(outvarnames[ivar]), VTKPointData()] = @view(qout[1:npoin,ivar])
        end

        # Optional caller-supplied nodal fields (see the NSD_2D writer).
        for (fname, fdata) in extra_fields
            vtkf[fname, VTKPointData()] = @view(fdata[1:npoin])
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
                      OUTPUT_DIR::String, inputs,
                      varnames, outvarnames,
                      outformat::HDF5;
                      nvar=1, qexact=zeros(1,nvar), case="",
                      μ_dsgs_pnode=nothing)
    
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
function write_output(SD, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs, varnames, outformat::HDF5; nvar=1, qexact=zeros(1,nvar), case="")
    
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
    
    MPI.Comm_rank(get_mpi_comm()) == 0 && println(string(" # Writing restart HDF5 file:", OUTPUT_DIR, "*.h5 ... DONE") )

end
function read_output(SD, INPUT_DIR::String, inputs, npoin, outformat::HDF5; nvar=1)
    
    #println(string(" # Reading restart HDF5 file:", INPUT_DIR, "*.h5 ...  ") )
    q, qe = read_hdf5(SD, INPUT_DIR, inputs, npoin, nvar)
    println_rank(string(" # Reading restart HDF5 file:", INPUT_DIR, "*.h5 ... DONE") ; msg_rank = rank)

    return q, qe
end


function write_hdf5(SD, mesh::St_mesh, q::AbstractArray, qe::AbstractArray, t, title::String, OUTPUT_DIR::String, inputs, varnames; iout=1, nvar=1, case="")
    # PERF: pull HDF5 into Jexpresso's namespace on first use; no-op
    # after that. Eager-loading HDF5 in src/Jexpresso.jl cost every
    # non-HDF5 run (city2d uses VTK output) tens of MB and seconds.
    _ensure_hdf5_loaded!()

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

function read_hdf5(SD, INPUT_DIR::String, inputs, npoin, nvar)
    _ensure_hdf5_loaded!()
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    mpi_size = MPI.Comm_size(comm)

    q  = zeros(Float64, npoin, nvar+1)
    qe = zeros(Float64, npoin, nvar+1)
    
    #read one HDF5 file time
    fout_name = string(INPUT_DIR, "/t.h5")
    time = rank == 0 ? convert(Float64, h5read(fout_name, "time")) : 0.0
    time = MPI.bcast(time, 0, comm)
    if inputs isa AbstractDict
        inputs[:tinit] = time
    else
        @warn "read_hdf5: cannot write :tinit into a NamedTuple inputs; " *
              "restart time = $time will be ignored unless you rebind " *
              "inputs in the caller."
    end
    #Write one HDF5 file per variable
    for ivar = 1:nvar
        fout_name   = string(INPUT_DIR, "/var_", ivar,"_",rank, ".h5")
        idx         = (ivar - 1)*npoin
        q[:, ivar]  = convert(Array{Float64, 1}, h5read(fout_name, "q"))
        qe[:, ivar] = convert(Array{Float64, 1}, h5read(fout_name, "qe"))
    end
    
    return q, qe
end

function write_NetCDF(SD::NSD_2D, mesh::St_mesh, q::Array, qaux::Array, mp,
                   connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original,
                   t, title::String, OUTPUT_DIR::String, inputs, varnames, outvarnames;
                   iout=1, nvar=1, qexact=zeros(1,nvar), case="")
    # PERF: pull NCDatasets into Jexpresso's namespace on first use.
    _ensure_netcdf_loaded!()

    if (isa(varnames, Tuple)    || isa(varnames, String) )   varnames    = collect(varnames) end
    if (isa(outvarnames, Tuple) || isa(outvarnames, String)) outvarnames = collect(outvarnames) end

    xx      = mesh.x
    yy      = mesh.y
    nvar    = size(varnames, 1)
    noutvar = max(nvar, size(outvarnames,1))

    ngr            = mesh.ngr
    ngl            = mesh.ngl
    npoin          = mesh.npoin
    nelem          = mesh.nelem
    nelem_semi_inf = mesh.nelem_semi_inf
    if (nelem_semi_inf > 0)
        subelem = Array{Int64}(undef, nelem*(ngl-1)^2+nelem_semi_inf*(ngl-1)*(ngr-1), 4)
    else
        subelem = Array{Int64}(undef, nelem*(ngl-1)^2, 4)
    end

    nsubelem = size(subelem, 1)
    isel = 1
    for iel = 1:nelem
        for i = 1:ngl-1
            for j = 1:ngl-1
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
    
    for iel = 1:nelem_semi_inf
        for i = 1:ngl-1
            for j = 1:ngr-1
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
    
    comm   = get_mpi_comm()
    rank   = MPI.Comm_rank(comm)
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
                      t, title::String, OUTPUT_DIR::String, inputs, varnames, outvarnames;
                      iout=1, nvar=1, qexact=zeros(1,nvar), case="")
    _ensure_netcdf_loaded!()

    if (isa(varnames, Tuple)    || isa(varnames, String) )   varnames    = collect(varnames) end
    if (isa(outvarnames, Tuple) || isa(outvarnames, String)) outvarnames = collect(outvarnames) end
    
    nvar    = size(varnames, 1)
    noutvar = max(nvar, size(outvarnames,1))
    npoin   = mesh.npoin
    nelem   = mesh.nelem
    ngl     = mesh.ngl
    
    nsubelem = mesh.nelem*(mesh.ngl-1)^3
    subelem  = Array{Int64}(undef, nsubelem, 8)
    
    isel = 1
    for iel = 1:nelem
        for i = 1:ngl-1
            for j = 1:ngl-1
                for k = 1:ngl-1
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
        nx[:]        = @view(mesh.coords[:,1])
        ny[:]        = @view(mesh.coords[:,2])
        nz[:]        = @view(mesh.coords[:,3])
        v2n[:]       = subelem
            
    end

end
