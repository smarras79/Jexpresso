using LinearSolve
using SnoopCompile
using WriteVTK
import SciMLBase

include("./plotting/jeplots.jl")

#----------------------------------------------------------------------------------------------------------------------------------------------
# ∂q/∂t = RHS -> q(x,t)
#----------------------------------------------------------------------------------------------------------------------------------------------
# PNG
function write_output(SD::NSD_1D, q::Matrix, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::PNG)

    #Reference values only (definied in initial conditions)
    
    nvar = length(varnames)
    for ivar = 1:nvar
        plot_results(SD, mesh, q[ivar,1:mesh.npoin], "initial", OUTPUT_DIR, varnames, inputs; iout=1, nvar=nvar, PT=nothing)
    end
end

function write_output(SD::NSD_1D, sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, varnames, outformat::PNG; nvar=1, qexact=zeros(1,nvar), case="")
    
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  "))
    
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
            plot_results!(SD, mesh, sol.u[iout][:], title, OUTPUT_DIR, varnames, inputs; iout=iout, nvar=nvar, fig=fig,color = color,p=p,marker=marker,PT=nothing)
            
        end
    else
        fig = Figure(size = (1200,800),fontsize=22)
        for iout = 1:size(sol.t[:], 1)
            title = string("sol.u at time ", sol.t[iout])
            plot_results(SD, mesh, sol.u[iout][:], title, OUTPUT_DIR, varnames, inputs; iout=iout, nvar=nvar,PT=nothing)
        end
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
            plot_triangulation(SD, mesh, sol.u[iout][:], title,  OUTPUT_DIR, inputs; iout=iout, nvar=nvar)
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

function write_output(sol::SciMLBase.LinearSolution, SD::NSD_2D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG; nvar=1)
    
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  ") )
    title = @sprintf "Solution to ∇⋅∇(q) = f"
    if inputs[:lplot_surf3d]
        plot_surf3d(SD, mesh, sol.u, title, OUTPUT_DIR; iout=1, nvar=1, smoothing_factor=inputs[:smoothing_factor])
    else
        plot_triangulation(SD, mesh, sol.u, title, OUTPUT_DIR, inputs;)
    end
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE") )
end


#------------
# HDF5 writer/reader
#------------

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


function write_hdf5(SD::NSD_2D, mesh::St_mesh, q::Array, qe::Array, title::String, OUTPUT_DIR::String, inputs::Dict, varnames; iout=1, nvar=1, case="")
    
    #Write one HDF5 file per variable
    for ivar = 1:nvar
        fout_name = string(OUTPUT_DIR, "/var_", ivar, ".h5")
        idx = (ivar - 1)*mesh.npoin
        h5write(fout_name, "q",  q[idx+1:ivar*mesh.npoin]);
        h5write(fout_name, "qe", qe[1:mesh.npoin, ivar]);

    end
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE") )
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
    #qe_temp = zeros(mesh.npoin,5)
    #qe_temp .= qexact
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
	q_new = zeros(new_size*4)
        q_exact1 = zeros(new_size,5)
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
                        mesh.connijk[e,l,m] = ip_new
                        mesh.poin_in_bdy_edge[iedge,k] = ip_new
		    elseif (unwind==1)
			#@info iter,ip,e,mesh.y[ip]
			ip_new = mesh.npoin + iter
                        mesh.x[ip_new] = xmax
                        mesh.connijk[e,l,m] = ip_new
                        mesh.y[ip_new] = mesh.y[ip]
			mesh.poin_in_bdy_edge[iedge,k] = ip_new
                        iter += 1
                	for ieq=1:4
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
			if(mesh.x[ip1] == xmax && mesh.y[ip1] == ymax)
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
    
    #npoin = mesh.npoin
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


    #Solution:
    fout_name = string(OUTPUT_DIR, "/iter_", iout, ".vtu")    
    vtkfile = vtk_grid(fout_name, mesh.x[1:npoin], mesh.y[1:npoin], mesh.y[1:npoin]*0.0, cells)
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
