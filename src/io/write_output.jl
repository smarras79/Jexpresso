using LinearSolve
using SnoopCompile
using WriteVTK
using HDF5

include("./plotting/jeplots.jl")

#----------------------------------------------------------------------------------------------------------------------------------------------
# ∂q/∂t = RHS -> q(x,t)
#----------------------------------------------------------------------------------------------------------------------------------------------
#
# PNG 1D
#
function write_output(SD::NSD_1D, q::Array, t, iout, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::PNG; nvar=1, qexact=zeros(1,nvar), case="")
    #OK
    nvar = length(varnames)
    qout = zeros(mesh.npoin)
    
    plot_results(SD, mesh, q[:], "initial", OUTPUT_DIR, varnames, inputs; iout=1, nvar=nvar, PT=nothing)
end

function write_output(SD::NSD_1D, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::PNG; nvar=1, qexact=zeros(1,nvar), case="")
    
    #println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  "))
    
    if (inputs[:plot_overlap])
        fig = Figure(size = (1200,800),fontsize=22)
        colors = ["Blue","Red","Green","Yellow","Black","Purple","Orange"]
        markers = [:circle, :rect, :diamond,:hexagon,:cross,:xcross,:utriangle,:dtriangle,:pentagon,:star4,:star8]
        p = []
        for iout = 1:size(sol.t[:], 1)
            icolor = mod(iout,size(colors,1))+1
            color = colors[icolor]
            imarker = mod(iout,size(markers,1))+1
            marker = markers[imarker]
            title = string("sol.u at time ", sol.t[iout])
            if (inputs[:backend] == CPU())
                plot_results!(SD, mesh, sol.u[iout][:], title, OUTPUT_DIR, varnames, inputs; iout=iout, nvar=nvar, fig=fig,color = color,p=p,marker=marker,PT=nothing)
            else
                uout = KernelAbstractions.allocate(CPU(),Float32, Int64(mesh.npoin))
                KernelAbstractions.copyto!(CPU(), uout, sol.u[iout][:])
                plot_results!(SD, mesh, uout, title, OUTPUT_DIR, varnames, inputs; iout=iout, nvar=nvar, fig=fig,color = color,p=p,marker=marker,PT=nothing)
            end
        end
    else
        fig = Figure(size = (1200,800),fontsize=22)
        for iout = 1:size(sol.t[:], 1)
            title = string("sol.u at time ", sol.t[iout])
            if (inputs[:backend] == CPU())
                plot_results(SD, mesh, sol.u[iout][:], title, OUTPUT_DIR, varnames, inputs; iout=iout, nvar=nvar,PT=nothing)
            else
                uout = KernelAbstractions.allocate(CPU(), TFloat, Int64(mesh.npoin*nvar))
                KernelAbstractions.copyto!(CPU(), uout, sol.u[iout][:])
                convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
                plot_results(SD, mesh, uout, title, OUTPUT_DIR, varnames, inputs; iout=iout, nvar=nvar,PT=nothing)
            end
        end
    end
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE ") )
end

#
# PNG 2D
#
function write_output(SD::NSD_2D, u::Array, t, iout, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::PNG; nvar=1, qexact=zeros(1,nvar), case="")

    #println(string(" # Writing 2D output to PNG file:", OUTPUT_DIR, "*.png ...  "))
    
    if inputs[:lplot_surf3d]
        for iout = 1:size(sol.t[:], 1)
            title = @sprintf "final solution at t=%6.4f" t
            plot_surf3d(SD, mesh, u[:], title, OUTPUT_DIR; iout=iout, nvar=nvar, smoothing_factor=inputs[:smoothing_factor])
        end
    else

        title = @sprintf "final solution at t=%6.4f" t
        if (inputs[:backend] == CPU())
            plot_triangulation(SD, mesh, u[:], title,  OUTPUT_DIR, inputs; iout=iout, nvar=nvar)
        else
            u = KernelAbstractions.allocate(CPU(), TFloat, Int64(mesh.npoin))
            KernelAbstractions.copyto!(CPU(),u, u[:])
            convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
            plot_triangulation(SD, mesh, u, title,  OUTPUT_DIR, inputs; iout=iout, nvar=nvar)
        end
    end
    println(string(" # Writing 2D output to PNG file:", OUTPUT_DIR, "*.png ...  DONE"))
end

#
function write_output(SD::NSD_2D, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::PNG; nvar=1, qexact=zeros(1,nvar), case="")

    #println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  "))
    
    if inputs[:lplot_surf3d]
        for iout = 1:size(sol.t[:], 1)
            title = @sprintf "final solution at t=%6.4f" sol.t[iout]
            plot_surf3d(SD, mesh, sol.u[iout][:], title, OUTPUT_DIR; iout=iout, nvar=nvar, smoothing_factor=inputs[:smoothing_factor])
        end
    else
        for iout = 1:size(sol.t[:],1)
            title = @sprintf "final solution at t=%6.4f" sol.t[iout]
            if (inputs[:backend] == CPU())
                plot_triangulation(SD, mesh, sol.u[iout][:], title,  OUTPUT_DIR, inputs; iout=iout, nvar=nvar)
            else
                u = KernelAbstractions.allocate(CPU(), TFloat, Int64(mesh.npoin))
                KernelAbstractions.copyto!(CPU(),u, sol.u[iout][:])
                convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
                plot_triangulation(SD, mesh, u, title,  OUTPUT_DIR, inputs; iout=iout, nvar=nvar)
            end
        end
    end
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE"))
end

#
# ASCII 2D
#
function write_output(SD::NSD_2D, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::ASCII; nvar=1, PT=nothing)
    
    #println(string(" # Writing output to ASCII file:", OUTPUT_DIR, "*.dat ...  ") )
    
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

#
# VTK 2D/3D
#
function write_output(SD, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::VTK; nvar=1, qexact=zeros(1,nvar), case="")
 
    #println(string(" # Writing output to VTK file:", OUTPUT_DIR, "*.vtu ...  ") )
    
    for iout = 1:size(sol.t[:],1)
        if (inputs[:backend] == CPU())
            title = @sprintf "final solution at t=%6.4f" sol.t[iout]
            write_vtk(SD, mesh, sol.u[iout][:], sol.t[iout], title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, qexact=qexact, case=case)
        else
            u = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin*nvar)
            KernelAbstractions.copyto!(CPU(),u,sol.u[iout][:])
            u_exact = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin,nvar+1)
            KernelAbstractions.copyto!(CPU(),u_exact,qexact)
            convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
            title = @sprintf "final solution at t=%6.4f" sol.t[iout]
            write_vtk(SD, mesh, u, sol.t[iout], title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, qexact=u_exact, case=case)
        end
    end
    println(string(" # Writing output to VTK file:", OUTPUT_DIR, "*.vtu ... DONE") )
end

function write_output(SD, u, t, iout, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::VTK; nvar=1, qexact=zeros(1,nvar), case="")
    
    title = @sprintf "final solution at t=%6.4f" iout
    if (inputs[:backend] == CPU())
        write_vtk(SD, mesh, u, t, title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, qexact=qexact, case=case)        
    else
        #VERIFY THIS on GPU
        u_gpu = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin*nvar)
        KernelAbstractions.copyto!(CPU(),u_gpu, u)
        u_exact = KernelAbstractions.allocate(CPU(),TFloat,mesh.npoin,nvar+1)
        KernelAbstractions.copyto!(CPU(),u_exact,qexact)
        convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
        write_vtk(SD, mesh, u_gpu, t, title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, qexact=u_exact, case=case)
        convert_mesh_arrays!(SD, mesh, inputs[:backend], inputs)

    end

    println(string(" # writing ", OUTPUT_DIR, "/iter", iout, ".vtu at t=", t, " s... DONE") )

end

# PNG 2D
#
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_2D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG; nvar=1)
    
    #println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  ") )
    
    title = @sprintf "Solution to ∇⋅∇(q) = f"
    if inputs[:lplot_surf3d]
        plot_surf3d(SD, mesh, sol.u, title, OUTPUT_DIR; iout=1, nvar=1, smoothing_factor=inputs[:smoothing_factor])
    else
        if (inputs[:backend] == CPU())
            plot_triangulation(SD, mesh, sol.u, title, OUTPUT_DIR, inputs;)
        else
            u = KernelAbstractions.allocate(CPU(), TFloat, Int64(mesh.npoin))
            KernelAbstractions.copyto!(CPU(),u, sol.u[:])
            convert_mesh_arrays_to_cpu!(SD, mesh, inputs)
            #@info u
            plot_triangulation(SD, mesh, u, title,  OUTPUT_DIR, inputs;)
        end
    end
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE") )
end


#------------
# HDF5 writer/reader
#------------
function write_output(SD, u::Array, t, iout, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::HDF5; nvar=1, qexact=zeros(1,nvar), case="")
    
    # println(string(" # Writing restart HDF5 file:", OUTPUT_DIR, "*.h5 ...  ") )
    iout = size(t,1)
    title = @sprintf "Final solution at t=%6.4f" t
    
    write_hdf5(SD, mesh, u, qexact, title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, case=case)
    
    println(string(" # Writing restart HDF5 file:", OUTPUT_DIR, "*.h5 ... DONE") )
    
end
function write_output(SD, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::HDF5; nvar=1, qexact=zeros(1,nvar), case="")
    
    #println(string(" # Writing restart HDF5 file:", OUTPUT_DIR, "*.h5 ...  ") )
    
    iout = size(sol.t[:],1)
    title = @sprintf "Final solution at t=%6.4f" sol.t[iout]

    write_hdf5(SD, mesh, sol.u[iout][:], qexact, title, OUTPUT_DIR, inputs, varnames; iout=iout, nvar=nvar, case=case)
    
    println(string(" # Writing restart HDF5 file:", OUTPUT_DIR, "*.h5 ... DONE") )
    
end
function read_output(SD::NSD_2D, INPUT_DIR::String, inputs::Dict, npoin, outformat::HDF5; nvar=1)
    
    #println(string(" # Reading restart HDF5 file:", INPUT_DIR, "*.h5 ...  ") )
    q, qe = read_hdf5(SD, INPUT_DIR, inputs, npoin, nvar)
    println(string(" # Reading restart HDF5 file:", INPUT_DIR, "*.h5 ... DONE") )

    return q, qe
end


function write_hdf5(SD, mesh::St_mesh, q::Array, qe::Array, title::String, OUTPUT_DIR::String, inputs::Dict, varnames; iout=1, nvar=1, case="")
    
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

#------------
# VTK writer
#------------
function write_vtk(SD::NSD_2D, mesh::St_mesh, q::Array, t, title::String, OUTPUT_DIR::String, inputs::Dict, varnames; iout=1, nvar=1, qexact=zeros(1,nvar), case="")

    outvars = varnames
    nvars = length(outvars)
    
    if (mesh.nelem_semi_inf > 0)
        subelem = Array{Int64}(undef, mesh.nelem*(mesh.ngl-1)^2+mesh.nelem_semi_inf*(mesh.ngl-1)*(mesh.ngr-1), 4)
        cells = [MeshCell(VTKCellTypes.VTK_QUAD, [1, 2, 4, 3]) for _ in 1:mesh.nelem*(mesh.ngl-1)^2+mesh.nelem_semi_inf*(mesh.ngl-1)*(mesh.ngr-1)]
    else
        subelem = Array{Int64}(undef, mesh.nelem*(mesh.ngl-1)^2, 4)
        cells = [MeshCell(VTKCellTypes.VTK_QUAD, [1, 2, 4, 3]) for _ in 1:mesh.nelem*(mesh.ngl-1)^2]
    end
    isel = 1
    npoin = mesh.npoin
    xx = zeros(size(mesh.x,1))
    yy = zeros(size(mesh.x,1))
    xx .= mesh.x
    yy .= mesh.y
    conn = zeros(mesh.nelem,mesh.ngl,mesh.ngl)
    conn .= mesh.connijk
    if ("Laguerre" in mesh.bdy_edge_type)
        conn_lag = zeros(mesh.nelem_semi_inf,mesh.ngl,mesh.ngr)
        conn_lag .= mesh.connijk_lag
    end
    poin_bdy = zeros(size(mesh.bdy_edge_type,1),mesh.ngl)
    poin_bdy .= mesh.poin_in_bdy_edge
    qe_temp = similar(qexact)
    
    if ("periodic1" in mesh.bdy_edge_type)
    	xmin = 1000000000.0
        ymax = -1000000000.0
        for e=1:mesh.nelem
            for i=1:mesh.ngl
                for j=1:mesh.ngl
                    ip = mesh.connijk[e,i,j]
                    ymax = max(ymax,mesh.y[ip])
                    xmin = min(xmin,mesh.x[ip])
                end
            end
        end
	xmax = -xmin
	nedges = size(mesh.bdy_edge_type,1)
        new_size = size(mesh.x,1)
        diff = new_size-mesh.npoin 
	q_new = zeros(new_size*nvar)
        q_exact1 = zeros(new_size,nvar+1)
        @info mesh.npoin, new_size, size(q_exact1), size(qexact)
        q_exact1[1:mesh.npoin,:] .= qexact[1:mesh.npoin,:] 
        
        for ieq = 1:nvars
	    ivar = new_size*(ieq-1)
            ivar1 = mesh.npoin*(ieq-1)
            #@info ivar,ivar1,new_size*ieq-diff,mesh.npoin*ieq
	    q_new[ivar+1:new_size*ieq-diff] .= q[ivar1+1:mesh.npoin*ieq]
	end
        iter = 1
    	for iedge = 1:nedges

    	    if (mesh.bdy_edge_type[iedge] == "periodic1")
		e = mesh.bdy_edge_in_elem[iedge]
		for k=1:mesh.ngl
                    ip = mesh.poin_in_bdy_edge[iedge,k]
		    xedge = mesh.x[ip]
		    unwind = 0
                    l = 0
		    m = 0
        	    dx = abs(mesh.x[mesh.connijk[e,2,1]]-mesh.x[mesh.connijk[e,mesh.ngl-1,1]])/(mesh.ngl-3)
		    for i=1:mesh.ngl
			for j=1:mesh.ngl
                            ip1 = mesh.connijk[e,i,j]
			    if (mesh.x[ip1] > xedge + (mesh.ngl)*dx)
				unwind=1
			    end
			    if (ip1 == ip)
				l=i
				m=j
			    end
			end
		    end
		    rep = 0
                    ip_rep = 0
        	    if (k == 1 || k == mesh.ngl)
			for ee=1:mesh.nelem
			    for i=1:mesh.ngl
				for j=1:mesh.ngl
				    ip1 = mesh.connijk[ee,i,j]
				    if (ip1 > mesh.npoin && mesh.y[ip1] == mesh.y[ip])
					ip_rep = ip1
					rep = 1
				    end
				end
			    end
			end
		    end
                    if (rep==1 && unwind==1)
			#@info iter,ip,e,mesh.y[ip] 
                        ip_new = ip_rep
                        if (l > 0 && m >0)
                            mesh.connijk[e,l,m] = ip_new
                        end
                        for iedge_1=1:nedges
                            if (iedge != iedge_1 && mesh.bdy_edge_in_elem[iedge_1] == e)
                                for i=1:mesh.ngl
                                    if (mesh.poin_in_bdy_edge[iedge_1,i] == ip)
                                            mesh.poin_in_bdy_edge[iedge_1,i] = ip_new
                                    end
                                end
                            end
                        end
                        mesh.poin_in_bdy_edge[iedge,k] = ip_new
		    elseif (unwind==1)
			#@info iter,ip,e,mesh.y[ip]
			ip_new = mesh.npoin + iter
                        mesh.x[ip_new] = xmax
                        if (l > 0 && m >0)
                            mesh.connijk[e,l,m] = ip_new
                        end
                        mesh.y[ip_new] = mesh.y[ip]
			for iedge_1=1:nedges
                            if (iedge != iedge_1 && mesh.bdy_edge_in_elem[iedge_1] == e)
                                for i=1:mesh.ngl
                                    if (mesh.poin_in_bdy_edge[iedge_1,i] == ip)
                                            mesh.poin_in_bdy_edge[iedge_1,i] = ip_new
                                    end
                                end
                            end
                        end
                        mesh.poin_in_bdy_edge[iedge,k] = ip_new
                        iter += 1
                	for ieq=1:nvar
        	            ivar = new_size*(ieq-1)
	                    ivar1 = mesh.npoin*(ieq-1)
                            q_new[ivar+ip_new] = q[ivar1+ip]
			end
			q_exact1[ip_new,:] .= qexact[ip,:]
		    end
		end
	    end
        end
        
        npoin += iter-1;
        if ("Laguerre" in mesh.bdy_edge_type)
	    
            e = mesh.nelem_semi_inf
	    iter = 1
	    for j=1:mesh.ngr
		ip = mesh.connijk_lag[e,mesh.ngl,j]
		if (j==1)
		    ip1=1
		    while (ip1 <= npoin)
                        #@info mesh.x[ip1], TFloat(xmax), mesh.y[ip1], TFloat(ymax)
                        if(mesh.x[ip1] == TFloat(xmax) && mesh.y[ip1] == TFloat(ymax))
			    ip_new = ip1
			end
			ip1 +=1
		    end
                    mesh.connijk_lag[e,mesh.ngl,j] = ip_new
		else
		    ip_new = npoin + iter
                    mesh.x[ip_new] = xmax
                    mesh.connijk_lag[e,mesh.ngl,j] = ip_new
                    mesh.y[ip_new] = mesh.y[ip]
                    iter += 1
		    for ieq=1:4
                        ivar = new_size*(ieq-1)
                        ivar1 = mesh.npoin*(ieq-1)
                        q_new[ivar+ip_new] = q[ivar1+ip]
                    end
		    q_exact1[ip_new,:] .= qexact[ip,:]
		end
	    end
	    npoin +=iter -1

	end
        q = q_new
	qexact = q_exact1
	
    end  

    if ("periodic2" in mesh.bdy_edge_type)
        xmax = 1000000000.0
        ymin = -1000000000.0
        for e=1:mesh.nelem
            for i=1:mesh.ngl
                for j=1:mesh.ngl
                    ip = mesh.connijk[e,i,j]
                    ymin = min(ymin,mesh.y[ip])
                    xmax = max(xmax,mesh.x[ip])
                end
            end
        end
        ymax = -ymin
        nedges = size(mesh.bdy_edge_type,1)
        new_size = size(mesh.x,1)
        diff = new_size-npoin
        q_new = zeros(new_size*nvar)
        q_exact1 = zeros(new_size,nvar+1)
        q_exact1[1:npoin,:] .= qexact[1:npoin,:]

        for ieq = 1:nvars
            ivar = new_size*(ieq-1)
            ivar1 = npoin*(ieq-1)
            #@info ivar,ivar1,new_size*ieq-diff,mesh.npoin*ieq
            q_new[ivar+1:new_size*ieq-diff] .= q[ivar1+1:npoin*ieq]
        end
        iter = 1
        for iedge = 1:nedges

            if (mesh.bdy_edge_type[iedge] == "periodic2")
                e = mesh.bdy_edge_in_elem[iedge]
                for k=1:mesh.ngl
                    ip = mesh.poin_in_bdy_edge[iedge,k]
                    yedge = mesh.y[ip]
                    unwind = 0
                    l = 0
                    m = 0
                    dy = abs(mesh.y[mesh.connijk[e,1,2]]-mesh.y[mesh.connijk[e,1,mesh.ngl-1]])/(mesh.ngl-3)
                    for i=1:mesh.ngl
                        for j=1:mesh.ngl
                            ip1 = mesh.connijk[e,i,j]
                            if (mesh.y[ip1] > yedge + (mesh.ngl)*dy)
                                unwind=1
                            end
                             if (ip1 == ip)
                                l=i
                                m=j
                            end
                        end
                    end
                    rep = 0
                    ip_rep = 0
                    if (k == 1 || k == mesh.ngl)
                        for ee=1:mesh.nelem
                            for i=1:mesh.ngl
                                for j=1:mesh.ngl
                                    ip1 = mesh.connijk[ee,i,j]
                                    if (ip1 > npoin && mesh.x[ip1] == mesh.x[ip])
                                        ip_rep = ip1
                                        rep = 1
                                    end
                                end
                            end
                        end
                    end
                    if (rep==1 && unwind==1)
                        #@info iter,ip,e,mesh.y[ip]

                        ip_new = ip_rep
                        if (l > 0 && m >0)
                            mesh.connijk[e,l,m] = ip_new
                        end
                        for iedge_1=1:nedges
                            if (iedge != iedge_1 && mesh.bdy_edge_in_elem[iedge_1] == e)
                                for i=1:mesh.ngl
                                    if (mesh.poin_in_bdy_edge[iedge_1,i] == ip)
                                            mesh.poin_in_bdy_edge[iedge_1,i] = ip_new
                                    end
                                end
                            end
                        end
                        mesh.poin_in_bdy_edge[iedge,k] = ip_new
                    elseif (unwind==1)
                        #@info iter,ip,e,mesh.y[ip]
                        ip_new = npoin + iter
                        mesh.y[ip_new] = ymax
                        if (l > 0 && m >0)
                            mesh.connijk[e,l,m] = ip_new
                        end
                        mesh.x[ip_new] = mesh.x[ip]
                        for iedge_1=1:nedges
                            if (iedge != iedge_1 && mesh.bdy_edge_in_elem[iedge_1] == e)
                                for i=1:mesh.ngl
                                    if (mesh.poin_in_bdy_edge[iedge_1,i] == ip)
                                            mesh.poin_in_bdy_edge[iedge_1,i] = ip_new
                                    end
                                end
                            end
                        end
                        mesh.poin_in_bdy_edge[iedge,k] = ip_new
                        iter += 1
                        for ieq=1:nvar
                            ivar = new_size*(ieq-1)
                            ivar1 = npoin*(ieq-1)
                            q_new[ivar+ip_new] = q[ivar1+ip]
                        end
                        q_exact1[ip_new,:] .= qexact[ip,:]
                    end
                end
            end
        end
        npoin += iter-1;
        q = q_new
        qexact = q_exact1

    end


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
    
    qout = copy(q)

    if (inputs[:CL] == CL())

        if (inputs[:SOL_VARS_TYPE] == TOTAL())
            
            #ρ
            qout[1:npoin] .= q[1:npoin]
            
            if (case == "rtb" || case == "mountain") && nvars >= 4
                
                #u = ρu/ρ
                ivar = 2
                idx = (ivar - 1)*npoin
                qout[idx+1:2*npoin] .= q[idx+1:2*npoin]./q[1:npoin]

                #v = ρv/ρ
                ivar = 3
                idx = (ivar - 1)*npoin
                qout[idx+1:3*npoin] .= q[idx+1:3*npoin]./q[1:npoin]
                if (size(qexact, 1) == npoin)

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

            for ivar = 2:nvars
                #u = ρu/ρ

                idx = (ivar - 1)*npoin
                qout[idx+1:ivar*npoin] .= (q[idx+1:ivar*npoin] .+ qexact[1:npoin,ivar])./(qout[1:npoin] .+ qexact[1:npoin,1]) .- qexact[1:npoin,ivar]./qexact[1:npoin,1]
                
                if (case == "rtb" || case == "mountain") && nvars >= 4
                    
                    if (size(qexact, 1) === npoin)
                        
                        ivar = 4
                        idx = (ivar - 1)*npoin
                        qout[idx+1:4*npoin] .= (q[idx+1:4*npoin] .+ qexact[1:npoin,4])./(qout[1:npoin] .+ qexact[1:npoin,1]) .- qexact[1:npoin,4]./qexact[1:npoin,1]
                    end
                end
            end
        end
        
    elseif (inputs[:CL] == NCL())
        
        #ρ
        qout[1:npoin] .= (q[1:npoin])

        for ivar = 1:nvars
            
            idx = (ivar - 1)*npoin
            qout[idx+1:ivar*npoin] .= (q[idx+1:ivar*npoin])

            if case === "rtb" && nvars >= 4
                
                if (inputs[:loutput_pert] == true && size(qexact, 1) === npoin)
                    outvars = ("dρ", "u", "v", "dθ")
                    #ρ'
                    qout[1:npoin] .= (q[1:npoin] .- qexact[1:npoin,1])
                    
                    #θ' = (ρθ - ρθref)/ρ = ρθ/ρ - ρrefθref/ρref
                    ivar = 4
                    idx = (ivar - 1)*npoin
                    qout[idx+1:4*npoin] .= (q[idx+1:4*npoin] .- qexact[1:npoin,4])
                    
                else
                    
                    #E = ρE/ρ
                    idx = 4*npoin
                    qout[idx+1:4*npoin] .= @views((q[2*npoin+1:4*npoin] .- 0.5*(q[npoin+1:2*npoin].*q[npoin+1:2*npoin] .+ q[npoin+1:3*npoin].*q[npoin+1:3*npoin])./q[1:npoin])./q[1:npoin]) #internal energy: p/((γ-1)ρ)
                end
            end
        end
    end
    
    if nvar > 4
        for ivar = 5:nvar
            idx = (ivar - 1)*npoin
            qout[idx+1:ivar*npoin] .= q[idx+1:ivar*npoin]
        end
    end

    #
    # Write solution:
    #
    fout_name = string(OUTPUT_DIR, "/iter_", iout, ".vtu")
    vtkfile = vtk_grid(fout_name, mesh.x[1:npoin], mesh.y[1:npoin], mesh.y[1:npoin]*TFloat(0.0), cells)
    for ivar = 1:nvar
        idx = (ivar - 1)*npoin
        vtkfile[string(varnames[ivar]), VTKPointData()] =  @view(qout[idx+1:ivar*npoin])
    end
    outfiles = vtk_save(vtkfile)
    mesh.x .= xx
    mesh.y .= yy
    mesh.connijk .= conn
    mesh.poin_in_bdy_edge .= poin_bdy
    qexact = copy(qe_temp)
    if ("Laguerre" in mesh.bdy_edge_type)
    	mesh.connijk_lag .= conn_lag 
    end
end

function write_vtk(SD::NSD_3D, mesh::St_mesh, q::Array, t, title::String, OUTPUT_DIR::String, inputs::Dict, varnames; iout=1, nvar=1, qexact=zeros(1,nvar), case="")
    outvars = varnames
    nvars = length(outvars)
    npoin = mesh.npoin
    xx = zeros(size(mesh.x,1))
    yy = zeros(size(mesh.x,1))
    zz = zeros(size(mesh.x,1))
    xx .= mesh.x 
    yy .= mesh.y
    zz .= mesh.z
    conn = zeros(mesh.nelem,mesh.ngl,mesh.ngl,mesh.ngl)
    conn .= mesh.connijk
    poin_bdy = zeros(size(mesh.bdy_face_type,1),mesh.ngl,mesh.ngl)
    poin_bdy .= mesh.poin_in_bdy_face
    qe_temp = similar(qexact)
    if ("periodic1" in mesh.bdy_face_type)
        xmax = mesh.xmax
        nfaces = size(mesh.bdy_face_type,1)
        new_size = size(mesh.x,1)
        diff = new_size-npoin
        q_new = zeros(new_size*nvar)
        q_exact1 = zeros(new_size,nvar+1)
        q_exact1[1:npoin,:] .= qexact[1:npoin,:]

        for ieq = 1:nvars
            ivar = new_size*(ieq-1)
            ivar1 = npoin*(ieq-1)
            q_new[ivar+1:new_size*ieq-diff] .= q[ivar1+1:npoin*ieq]
        end
        iter = 1
        for iface = 1:nfaces

            if (mesh.bdy_face_type[iface] == "periodic1")
                e = mesh.bdy_face_in_elem[iface]
                for k=1:mesh.ngl
                    for l=1:mesh.ngl
                        ip = mesh.poin_in_bdy_face[iface,k,l]
                        xface = mesh.x[ip]
                        unwind = 0
                        ll = 0
                        m = 0
                        n = 0
                        dx = abs(mesh.x[mesh.connijk[e,2,1,1]]-mesh.x[mesh.connijk[e,mesh.ngl-1,1,1]])/(mesh.ngl-3)
                        for i=1:mesh.ngl
                            for j=1:mesh.ngl
                                for kk=1:mesh.ngl
                                    ip1 = mesh.connijk[e,i,j,kk]
                                    if (mesh.x[ip1] > xface + (mesh.ngl)*dx)
                                        unwind=1
                                    end
                                    if (ip1 == ip)
                                        ll=i
                                        m=j
                                        n = kk
                                    end
                                end
                            end
                        end
                        rep = 0
                        ip_rep = 0
                        if (k == 1 || k == mesh.ngl || l == 1 || l == mesh.ngl)
                            for ee=1:mesh.nelem
                                for i=1:mesh.ngl
                                    for j=1:mesh.ngl
                                        for kk=1:mesh.ngl
                                            ip1 = mesh.connijk[ee,i,j,kk]
                                            if (ip1 > npoin && mesh.y[ip1] == mesh.y[ip] && mesh.z[ip1] == mesh.z[ip])
                                                ip_rep = ip1
                                                rep = 1
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        if (rep==1 && unwind==1)
                            ip_new = ip_rep
                            if (ll > 0 && m > 0 && n > 0)
                                mesh.connijk[e,ll,m,n] = ip_new
                            end
                            for iface_1 = 1:nfaces
                                if (iface != iface_1 && mesh.bdy_face_in_elem[iface_1] == e)
                                    for i=1:mesh.ngl
                                        for j =1:mesh.ngl
                                            if (mesh.poin_in_bdy_face[iface_1,i,j] == ip)
                                                mesh.poin_in_bdy_face[iface_1,i,j] = ip_new
                                            end
                                        end
                                    end
                                end
                            end
                            mesh.poin_in_bdy_face[iface,k,l] = ip_new
                        elseif (unwind==1)
                            ip_new = npoin + iter
                            mesh.x[ip_new] = xmax
                            if (ll > 0 && m > 0 && n > 0)
                                mesh.connijk[e,ll,m,n] = ip_new
                            end
                            mesh.y[ip_new] = mesh.y[ip]
                            mesh.z[ip_new] = mesh.z[ip]
                            for iface_1 = 1:nfaces
                                if (iface != iface_1 && mesh.bdy_face_in_elem[iface_1] == e)
                                    for i=1:mesh.ngl
                                        for j =1:mesh.ngl
                                            if (mesh.poin_in_bdy_face[iface_1,i,j] == ip)
                                                mesh.poin_in_bdy_face[iface_1,i,j] = ip_new
                                            end
                                        end 
                                    end 
                                end 
                            end
                            mesh.poin_in_bdy_face[iface,k,l] = ip_new
                            iter += 1
                            for ieq=1:nvar
                                ivar = new_size*(ieq-1)
                                ivar1 = npoin*(ieq-1)
                                q_new[ivar+ip_new] = q[ivar1+ip]
                            end
                            q_exact1[ip_new,:] .= qexact[ip,:]
                        end
                    end
                end
            end
        end
        npoin += iter-1;
        q = q_new
        qexact = q_exact1
    end
    
    if ("periodic2" in mesh.bdy_face_type)
        zmax = mesh.zmax
        nfaces = size(mesh.bdy_face_type,1)
        new_size = size(mesh.x,1)
        diff = new_size-npoin
        q_new = zeros(new_size*nvar)
        q_exact1 = zeros(new_size,nvar+1)
        q_exact1[1:npoin,:] .= qexact[1:npoin,:]

        for ieq = 1:nvars
            ivar = new_size*(ieq-1)
            ivar1 = npoin*(ieq-1)
            q_new[ivar+1:new_size*ieq-diff] .= q[ivar1+1:npoin*ieq]
        end
        iter = 1
        for iface = 1:nfaces

            if (mesh.bdy_face_type[iface] == "periodic2")
                e = mesh.bdy_face_in_elem[iface]
                for k=1:mesh.ngl
                    for l=1:mesh.ngl
                        ip = mesh.poin_in_bdy_face[iface,k,l]
                        zface = mesh.z[ip]
                        unwind = 0
                        ll = 0
                        m = 0
                        n = 0
                        dz = abs(mesh.z[mesh.connijk[e,1,1,2]]-mesh.z[mesh.connijk[e,1,1,mesh.ngl-1]])/(mesh.ngl-3)
                        for i=1:mesh.ngl
                            for j=1:mesh.ngl
                                for kk=1:mesh.ngl
                                    ip1 = mesh.connijk[e,i,j,kk]
                                    if (mesh.z[ip1] > zface + (mesh.ngl)*dz)
                                        unwind=1
                                    end
                                    if (ip1 == ip)
                                        ll=i
                                        m=j
                                        n = kk
                                    end
                                end
                            end
                        end
                        rep = 0
                        ip_rep = 0
                        if (k == 1 || k == mesh.ngl || l == 1 || l == mesh.ngl)
                            for ee=1:mesh.nelem
                                for i=1:mesh.ngl
                                    for j=1:mesh.ngl
                                        for kk=1:mesh.ngl
                                            ip1 = mesh.connijk[ee,i,j,kk]
                                            if (ip1 > npoin && mesh.x[ip1] == mesh.x[ip] && mesh.y[ip1] == mesh.y[ip])
                                                ip_rep = ip1
                                                rep = 1
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        if (rep==1 && unwind==1)
                            ip_new = ip_rep
                            if (ll > 0 && m > 0 && n > 0)
                                mesh.connijk[e,ll,m,n] = ip_new
                            end
                            for iface_1 = 1:nfaces
                                if (iface != iface_1 && mesh.bdy_face_in_elem[iface_1] == e)
                                    for i=1:mesh.ngl
                                        for j =1:mesh.ngl
                                            if (mesh.poin_in_bdy_face[iface_1,i,j] == ip)
                                                mesh.poin_in_bdy_face[iface_1,i,j] = ip_new
                                            end
                                        end
                                    end
                                end
                            end
                            mesh.poin_in_bdy_face[iface,k,l] = ip_new
                        elseif (unwind==1)
                            #@info iter,ip,e,mesh.y[ip]
                            ip_new = npoin + iter
                            mesh.z[ip_new] = zmax
                            if (ll > 0 && m > 0 && n > 0)
                                mesh.connijk[e,ll,m,n] = ip_new
                            end
                            mesh.y[ip_new] = mesh.y[ip]
                            mesh.x[ip_new] = mesh.x[ip]
                            for iface_1 = 1:nfaces
                                if (iface != iface_1 && mesh.bdy_face_in_elem[iface_1] == e)
                                    for i=1:mesh.ngl
                                        for j =1:mesh.ngl
                                            if (mesh.poin_in_bdy_face[iface_1,i,j] == ip)
                                                mesh.poin_in_bdy_face[iface_1,i,j] = ip_new
                                            end
                                        end 
                                    end 
                                end 
                            end
                            mesh.poin_in_bdy_face[iface,k,l] = ip_new
                            iter += 1
                            for ieq=1:nvar
                                ivar = new_size*(ieq-1)
                                ivar1 = npoin*(ieq-1)
                                q_new[ivar+ip_new] = q[ivar1+ip]
                            end
                            q_exact1[ip_new,:] .= qexact[ip,:]
                        end
                    end
                end
            end
        end
        npoin += iter-1;
        q = q_new
        qexact = q_exact1
    end
    
    if ("periodic3" in mesh.bdy_face_type)
        ymax = mesh.ymax
        nfaces = size(mesh.bdy_face_type,1)
        new_size = size(mesh.x,1)
        diff = new_size-npoin
        q_new = zeros(new_size*nvar)
        q_exact1 = zeros(new_size,nvar+1)
        q_exact1[1:npoin,:] .= qexact[1:npoin,:]

        for ieq = 1:nvars
            ivar = new_size*(ieq-1)
            ivar1 = npoin*(ieq-1)
            q_new[ivar+1:new_size*ieq-diff] .= q[ivar1+1:npoin*ieq]
        end
        iter = 1
        for iface = 1:nfaces

            if (mesh.bdy_face_type[iface] == "periodic3")
                e = mesh.bdy_face_in_elem[iface]
                for k=1:mesh.ngl
                    for l=1:mesh.ngl
                        ip = mesh.poin_in_bdy_face[iface,k,l]
                        yface = mesh.y[ip]
                        unwind = 0
                        ll = 0
                        m = 0
                        n = 0
                        dy = abs(mesh.y[mesh.connijk[e,1,2,1]]-mesh.y[mesh.connijk[e,1,mesh.ngl-1,1]])/(mesh.ngl-3)
                        for i=1:mesh.ngl
                            for j=1:mesh.ngl
                                for kk=1:mesh.ngl
                                    ip1 = mesh.connijk[e,i,j,kk]
                                    if (mesh.y[ip1] > yface + (mesh.ngl)*dy)
                                        unwind=1
                                    end
                                    if (ip1 == ip)
                                        ll=i
                                        m=j
                                        n=kk
                                    end
                                end
                            end
                        end
                        rep = 0
                        ip_rep = 0
                        if (k == 1 || k == mesh.ngl || l == 1 || l == mesh.ngl)
                            for ee=1:mesh.nelem
                                for i=1:mesh.ngl
                                    for j=1:mesh.ngl
                                        for kk=1:mesh.ngl
                                            ip1 = mesh.connijk[ee,i,j,kk]
                                            if (ip1 > npoin && mesh.x[ip1] == mesh.x[ip] && mesh.z[ip1] == mesh.z[ip])
                                                ip_rep = ip1
                                                rep = 1
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        if (rep==1 && unwind==1)
                            ip_new = ip_rep
                            if (ll > 0 && m > 0 && n > 0)
                                mesh.connijk[e,ll,m,n] = ip_new
                            end
                            for iface_1 = 1:nfaces
                                if (iface != iface_1 && mesh.bdy_face_in_elem[iface_1] == e)
                                    for i=1:mesh.ngl
                                        for j =1:mesh.ngl
                                            if (mesh.poin_in_bdy_face[iface_1,i,j] == ip)
                                                mesh.poin_in_bdy_face[iface_1,i,j] = ip_new
                                            end
                                        end
                                    end
                                end
                            end
                            mesh.poin_in_bdy_face[iface,k,l] = ip_new
                        elseif (unwind==1)
                            ip_new = npoin + iter
                            mesh.y[ip_new] = ymax
                            if (ll > 0 && m > 0 && n > 0)
                                mesh.connijk[e,ll,m,n] = ip_new
                            end

                            mesh.z[ip_new] = mesh.z[ip]
                            mesh.x[ip_new] = mesh.x[ip]
                            for iface_1 = 1:nfaces
                                if (iface != iface_1 && mesh.bdy_face_in_elem[iface_1] == e)
                                    for i=1:mesh.ngl
                                        for j =1:mesh.ngl
                                            if (mesh.poin_in_bdy_face[iface_1,i,j] == ip)
                                                mesh.poin_in_bdy_face[iface_1,i,j] = ip_new
                                            end
                                        end 
                                    end 
                                end 
                            end
                            mesh.poin_in_bdy_face[iface,k,l] = ip_new
                            iter += 1
                            for ieq=1:nvar
                                ivar = new_size*(ieq-1)
                                ivar1 = npoin*(ieq-1)
                                q_new[ivar+ip_new] = q[ivar1+ip]
                            end
                            q_exact1[ip_new,:] .= qexact[ip,:]
                        end
                    end
                end
            end
        end
        npoin += iter-1;
        q = q_new
        qexact = q_exact1
    end

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

    npoin1 = mesh.npoin
    mesh.npoin = npoin
    #npoin = mesh.npoin
    qout = copy(q)
    TF = eltype(q)
    ρ = zeros(TF, mesh.npoin,1)
    u = zeros(TF, mesh.npoin,1)
    v = zeros(TF, mesh.npoin,1)
    w = zeros(TF, mesh.npoin,1)
    e_tot = zeros(TF, mesh.npoin,1)
    qt = zeros(TF, mesh.npoin,1)
    ql = zeros(TF, mesh.npoin,1)
    θ = zeros(TF, mesh.npoin,1)
    param_set = create_updated_TD_Parameters(TF(101325.0))
    _grav = TF(TP.grav(param_set))


    if (inputs[:CL] == CL())
        if (inputs[:SOL_VARS_TYPE] == TOTAL())
            #ρ
            qout[1:mesh.npoin] .= q[1:mesh.npoin]
            ρ[:] .= q[1:mesh.npoin]
            
            if (case == "rtb" || case == "mountain") && nvars >= 4
                #u = ρu/ρ
                ivar = 2
                idx = (ivar - 1)*mesh.npoin
                qout[idx+1:ivar*mesh.npoin] .= q[idx+1:ivar*mesh.npoin]./q[1:mesh.npoin]
                u[:] .= qout[idx+1:ivar*mesh.npoin]

                #v = ρv/ρ
                ivar = 3
                idx = (ivar - 1)*mesh.npoin
                qout[idx+1:ivar*mesh.npoin] .= q[idx+1:ivar*mesh.npoin]./q[1:mesh.npoin]
                v[:] .= qout[idx+1:ivar*mesh.npoin]
                
                #w = ρw/ρ
                ivar = 4
                idx = (ivar - 1)*mesh.npoin
                qout[idx+1:ivar*mesh.npoin] .= q[idx+1:ivar*mesh.npoin]./q[1:mesh.npoin]
                w[:] .= qout[idx+1:ivar*mesh.npoin]
                
                if (size(qexact, 1) == npoin)

                    if inputs[:loutput_pert] == true
                        #ρ'
                        qout[1:npoin] .= q[1:npoin] .- qexact[1:npoin,1]
                        
                        #θ' = (ρθ - ρθref)/ρ = ρθ/ρ - ρrefθref/ρref
                        ivar = 5
                        idx = (ivar - 1)*npoin
                        qout[idx+1:5*npoin] .= q[idx+1:5*npoin]./q[1:npoin] .- qexact[1:npoin,5]./qexact[1:npoin,1]
                    else
                        ivar = 5
                        idx = (ivar - 1)*mesh.npoin
                        qout[idx+1:5*mesh.npoin] .= q[idx+1:5*mesh.npoin]./q[1:mesh.npoin]
                        e_tot[:] .= qout[idx+1:ivar*mesh.npoin]

                    end
                end
                
                if (inputs[:lbomex])
                    #qt = ρqt/ρ
                    ivar = 6
                    idx = (ivar - 1)*mesh.npoin
                    qout[idx+1:ivar*mesh.npoin] .= q[idx+1:ivar*mesh.npoin]./q[1:mesh.npoin]
                    qt[:] .= qout[idx+1:ivar*mesh.npoin]


                    #ql = ρql/ρ
                    ivar = 7
                    idx = (ivar - 1)*mesh.npoin
                    qout[idx+1:ivar*mesh.npoin] .= q[idx+1:ivar*mesh.npoin]./q[1:mesh.npoin]
                    ql[:] .= qout[idx+1:ivar*mesh.npoin]
                    
                    for i = 1:mesh.npoin
                        e_kin = 0.5 * (u[i]^2 + v[i]^2 + w[i]^2)
                        e_pot = _grav * mesh.z[i]
                        e_int::TF = e_tot[i] - e_pot - e_kin
                        q_pt = TD.PhasePartition(qt[i], ql[i])
                        Temp = TD.air_temperature(param_set, e_int, q_pt)
                        θ[i] = TD.virtual_pottemp(param_set, Temp, ρ[i], q_pt)
                    end
                end
            end
        else
            qout[1:npoin] .= q[1:npoin]

            for ivar = 2:nvars
                #u = ρu/ρ

                idx = (ivar - 1)*npoin
                qout[idx+1:ivar*npoin] .= (q[idx+1:ivar*npoin] .+ qexact[1:npoin,ivar])./(qout[1:npoin] .+ qexact[1:npoin,1]) .- qexact[1:npoin,ivar]./qexact[1:npoin,1]

                if (case == "rtb" || case == "mountain") && nvars >= 4

                    if (size(qexact, 1) === npoin)

                        ivar = 5
                        idx = (ivar - 1)*npoin
                        qout[idx+1:5*npoin] .= (q[idx+1:5*npoin] .+ qexact[1:npoin,5])./(qout[1:npoin] .+ qexact[1:npoin,1]) .- qexact[1:npoin,5]./qexact[1:npoin,1]
                    end
                end
            end

        end
    end

    #Solution:
    fout_name = string(OUTPUT_DIR, "/iter_", iout, ".vtu")
    vtkfile = vtk_grid(fout_name, mesh.x[1:mesh.npoin], mesh.y[1:mesh.npoin], mesh.z[1:mesh.npoin], cells)
    for ivar = 1:nvars
        idx = (ivar - 1)*npoin
        vtkfile[string(varnames[ivar]), VTKPointData()] =  @view(qout[idx+1:ivar*npoin])
    end
    if (inputs[:lbomex])
        vtkfile["theta", VTKPointData()] =  @view(θ[:])
    end
    outfiles = vtk_save(vtkfile)
    mesh.npoin = npoin1
    mesh.x .= xx
    mesh.y .= yy
    mesh.z .= zz
    mesh.connijk .= conn
    mesh.poin_in_bdy_face .= poin_bdy
    qexact = copy(qe_temp)
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
end
