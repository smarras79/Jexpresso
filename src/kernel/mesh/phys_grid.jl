Base.@kwdef mutable struct phys_grid{T <: AbstractFloat, dims1, dims2, dims3, dims4, dims5, dims6, dims7, backend}

    x = KernelAbstractions.zeros(backend,T, dims1)
    y = KernelAbstractions.zeros(backend,T, dims2)
    x_col = KernelAbstractions.zeros(backend,T, dims3)
    y_col = KernelAbstractions.zeros(backend,T, dims3)
    z = KernelAbstractions.zeros(backend,T, dims4)

    ncol = dims3[1]
    nlev = dims4[1]
    nx = dims1[1]
    ny = dims2[1]

    p     = KernelAbstractions.zeros(backend,T, dims4)
    t     = KernelAbstractions.zeros(backend,T, dims4)
    p_lay = KernelAbstractions.zeros(backend,T, dims5)
    t_lay = KernelAbstractions.zeros(backend,T, dims5)
    qv    = KernelAbstractions.zeros(backend,T, dims4)
    qc    = KernelAbstractions.zeros(backend,T, dims4)
    qi    = KernelAbstractions.zeros(backend,T, dims4)
    ρ     = KernelAbstractions.zeros(backend,T, dims4)
    qv_lay    = KernelAbstractions.zeros(backend,T, dims5)
    qc_lay    = KernelAbstractions.zeros(backend,T, dims5)
    qi_lay    = KernelAbstractions.zeros(backend,T, dims5)
    ρ_lay     = KernelAbstractions.zeros(backend,T, dims5)

    dry   = KernelAbstractions.zeros(backend,T, dims5)
    rel   = KernelAbstractions.zeros(backend,T, dims5)
    vmr   = KernelAbstractions.zeros(backend,T, dims4)
    vmr_lay = KernelAbstractions.zeros(backend,T, dims5)
    
    el_correspondance = KernelAbstractions.zeros(backend, TInt, dims4)
    ω_dyn = KernelAbstractions.zeros(backend, T, dims7)
    int_points = KernelAbstractions.zeros(backend, TInt, dims6)
    ω_x = KernelAbstractions.zeros(backend,T, dims1)
    ω_y = KernelAbstractions.zeros(backend,T, dims2)
    ω_z = KernelAbstractions.zeros(backend,T, dims4)
end

function init_phys_grid(mesh,inputs,nlay,nx,ny,xmin,xmax,ymin,ymax,zmin,zmax,backend)
    
    ncol = nx*ny
    dims1 = (Int64(nx))
    dims2 = (Int64(ny))
    dims3 = (Int64(ncol))
    dims4 = (Int64(nlay+1),Int64(ncol))
    dims5 = (Int64(nlay), Int64(ncol))
    dims6 = (Int64(nlay+1),Int64(ncol), 8)
    dims7 = (Int64(nlay+1), Int64(ncol), Int64(mesh.ngl), Int64(mesh.ngl), Int64(mesh.ngl))
    Rad_grid = phys_grid{TFloat, dims1, dims2, dims3, dims4, dims5, dims6, dims7, backend}()
    dx = (xmax - xmin)/nx
    dy = (ymax - ymin)/ny
    dz = (zmax - zmin)/(nlay+1)
    ξωx  = basis_structs_ξ_ω!(inputs[:interpolation_nodes], nx-1, inputs[:backend])
    ξx = ξωx.ξ
    ξωy  = basis_structs_ξ_ω!(inputs[:interpolation_nodes], ny-1, inputs[:backend])
    ξy = ξωy.ξ
    ξωz  = basis_structs_ξ_ω!(inputs[:interpolation_nodes], nlay, inputs[:backend])
    ξz = ξωz.ξ
    scalex = (xmax - xmin)/2
    scaley = (ymax - ymin)/2
    scalez = (zmax - zmin)/2
    if (backend == CPU())
        for ix=1:nx
            Rad_grid.x[ix] = scalex*ξx[ix]#xmin + (ix-1)*dx
            for iy=1:ny
                Rad_grid.y[iy] = scaley*ξy[iy]#ymin + (iy-1)*dy
                icol = ix + nx*(iy-1)
                Rad_grid.x_col[icol] = scalex*ξx[ix]#xmin + (ix-1)*dx
                Rad_grid.y_col[icol] = scaley*ξy[iy]#ymin + (iy-1)*dy
            end
        end
        for icol = 1:ncol
            for ilev = 1:nlay+1
                Rad_grid.z[ilev,icol] = scalez*(ξz[ilev]+1)#zmin + (ilev-1)*dz
            end
        end
        if (inputs[:lwarp])
            warp_phys_grid!(Rad_grid.x_col,Rad_grid.y_col,Rad_grid.z,ncol,nlay)
        end
        for icol=1:ncol
            Rad_grid.ω_z[:,icol] = BarycentricWeights(Rad_grid.z[:,icol])
        end
        Rad_grid.ω_x = BarycentricWeights(Rad_grid.x)
        Rad_grid.ω_y = BarycentricWeights(Rad_grid.y)
        store_element_maxima!(mesh,mesh.el_max,mesh.el_min,mesh.nelem,mesh.ngl)
        store_mesh_to_phys_grid_correspondance!(Rad_grid,mesh.el_max,mesh.el_min,mesh.nelem,Rad_grid.ncol,Rad_grid.nlev-1,mesh)
    else
        k = build_phys_grid_gpu!(backend,Int64(ncol),Int64(nlay+1))
        k(Rad_grid.x, Rad_grid.y, Rad_grid.z, nx, ny, dx, dy, dz, xmin, ymin, zmin;
          ndrange = ncol*(nlay+1), workgroupsize= nlay+1)
        x_col = KernelAbstractions.allocate(CPU(), TFloat, size(Rad_grid.x_col))
        KernelAbstractions.copyto!(CPU(),x_col,Rad_grid.x_col)
        y_col = KernelAbstractions.allocate(CPU(), TFloat, size(Rad_grid.y_col))
        KernelAbstractions.copyto!(CPU(),y_col,Rad_grid.y_col)
        z = KernelAbstractions.allocate(CPU(), TFloat, size(Rad_grid.z))
        KernelAbstractions.copyto!(CPU(),z,Rad_grid.z)
        warp_phys_grid!(x_col,y_col,z,ncol,nlay)
        KernelAbstractions.copyto!(backend, Rad_grid.z, z)
        k = BarycentricWeights_gpu!(backend,Int64(nx))
        k(Rad_grid.x, nx, Rad_grid.ω_x; ndrange = nx)
        k = BarycentricWeights_gpu!(backend,Int64(ny))
        k(Rad_grid.y, ny, Rad_grid.ω_y; ndrange = ny)
        for icol=1:ncol
            k = BarycentricWeights_gpu!(backend,Int64(nlay+1))
            k(Rad_grid.z[:,icol], nlay+1, Rad_grid.ω_z[:,icol]; ndrange = nlay+1)
        end
    end
    return Rad_grid
end

@kernel function build_phys_grid_gpu!(x,y,z,nx,ny,dx,dy,dz,xmin,ymin,zmin)
    
    icol = @index(Group,Linear)
    ilev = @index(Local,Linear)
    T = eltype(x)
    ix = T(mod(icol,nx)) 
    iy = T((icol - ix)/nx + T(1))
    x[icol] = xmin + (ix-1)*dx
    y[icol] = ymin + (iy-1)*dy
    z[ilev,icol] = zmin + (ilev-1)*dz
    
end

function store_element_maxima!(mesh,el_max,el_min,nelem,ngl)

    for e = 1:nelem
        xmin = 100000000
        xmax = -100000000
        ymin = 100000000
        ymax = -100000000
        zmin = 100000000
        zmax = -100000000
        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip = mesh.connijk[e,i,j,k]
                    xmin = min(xmin,mesh.x[ip])
                    xmax = max(xmax,mesh.x[ip])
                    ymin = min(ymin,mesh.y[ip])
                    ymax = max(ymax,mesh.y[ip])
                    zmin = min(zmin,mesh.z[ip])
                    zmax = max(zmax,mesh.z[ip])
                end
            end
        end
        el_max[e,1] = xmax
        el_max[e,2] = ymax
        el_max[e,3] = zmax
        el_min[e,1] = xmin
        el_min[e,2] = ymin
        el_min[e,3] = zmin
    end
end

@kernel function  store_element_maxima_gpu!(x,y,z,connijk,el_max,el_min)
    e = @index(Group,Linear)
    il = @index(Local,NTuple)
    @inbounds i_x = il[1]
    @inbounds i_y = il[2]
    @inbounds i_z = il[3]
    ip = connijk[e,i_x,i_y,i_z]
    @inbounds KernelAbstractions.@atomic el_min[e,1] = min(el_min[e,1],mesh.x[ip])
    @inbounds KernelAbstractions.@atomic el_max[e,1] = max(el_max[e,1],mesh.x[ip])
    @inbounds KernelAbstractions.@atomic el_min[e,2] = min(el_min[e,2],mesh.y[ip])
    @inbounds KernelAbstractions.@atomic el_max[e,2] = max(el_max[e,2],mesh.y[ip])
    @inbounds KernelAbstractions.@atomic el_min[e,3] = min(el_min[e,3],mesh.x[ip])
    @inbounds KernelAbstractions.@atomic el_max[e,3] = max(el_max[e,3],mesh.z[ip])
end

function store_mesh_to_phys_grid_correspondance!(phys_grid,el_max,el_min,nelem,ncol,nlay,mesh)

    for ilay=1:nlay+1
        for icol=1:ncol
            x = phys_grid.x_col[icol]
            y = phys_grid.y_col[icol]
            z = phys_grid.z[ilay,icol]
            e = 1
            found = false
            while (e <= nelem && found == false)
                if (x >= el_min[e,1] && y >= el_min[e,2] && z >= el_min[e,3] && x <= el_max[e,1] && y <= el_max[e,2] && z <= el_max[e,3])
                    phys_grid.el_correspondance[ilay,icol] = e
                    found = true
                end
                e += 1
            end
            if (found == false)
                @info "error could not find corresponding element on dynamics mesh, something is broken"
                @info e, x, y, z, maximum(el_max[:,1]), minimum(el_min[:,1])
            end
            S = 0
            e = phys_grid.el_correspondance[ilay,icol]
            match_found = false
            for i=1:mesh.ngl
                for j=1:mesh.ngl
                    for k=1:mesh.ngl
                        ip = mesh.connijk[e,i,j,k]
                        if (AlmostEqual(x, mesh.x[ip]) && AlmostEqual(y, mesh.y[ip]) && AlmostEqual(z, mesh.z[ip]))
                            match_found = true
                            phys_grid.ω_dyn[ilay,icol,i,j,k] = 1
                        else 
                            phys_grid.ω_dyn[ilay,icol,i,j,k] = 1/ (sqrt((mesh.x[ip] - x)^2 + (mesh.y[ip] - y)^2 + (mesh.z[ip] - z)^2))
                        end
                    end
                end
            end
            if (match_found)
                for i=1:mesh.ngl
                    for j=1:mesh.ngl
                        for k=1:mesh.ngl
                            if (phys_grid.ω_dyn[ilay,icol,i,j,k] != 1)
                                phys_grid.ω_dyn[ilay,icol,i,j,k] = 0
                            end
                        end
                    end
                end
            end
            S = sum(phys_grid.ω_dyn[ilay,icol,:,:,:])
            phys_grid.ω_dyn[ilay,icol,:,:,:] .= phys_grid.ω_dyn[ilay,icol,:,:,:]/S
        end
    end
    ##storing weight before periodicity is applied
end

@kernel function store_mesh_to_phys_grid_correspondance_gpu!(el_correspondance,el_max,el_min,nelem,ncol,nlay)
    icol = @index(Group,Linear)
    ilev = @index(Local,Linear)
    T = eltype(nelem)
    e = T(1)
    found = false
    while (e <= nelem && found == false)
        if (x> el_min[e,1] && y > el_min[e,2] && z > el_min[e,3] && x < el_max[e,1] && y < el_max[e,2] && z < el_max[e,3])
            el_correspondance[ilev,icol] = e
            found = true
        else
            e+=1
        end
    end

end

function interpolate_to_phys_grid!(mesh,phys_grid,uaux,qe,nlay,ncol,P,T,qc,qi,ρ,lpert)

    PhysConst = PhysicalConst{Float64}()
    A = KernelAbstractions.zeros(CPU(), TFloat, mesh.npoin,6)
    ω = KernelAbstractions.zeros(CPU(), TFloat,mesh.ngl,mesh.ngl,mesh.ngl)
    B = KernelAbstractions.zeros(CPU(), TFloat, 6)
    A[:,2] .= uaux[:,end]
    A[:,3] .= T
    A[:,4] .= qc
    A[:,5] .= qi
    if (lpert)
        A[:,1] .= uaux[:,1] .+ qe[:,1]
        A[:,6] .= (uaux[:,6] .+ qe[:,6])./(A[:,1]) .- qc .- qi
    else
        A[:,1] .= uaux[:,1]
        A[:,6] .= uaux[:,6] .+ qe[:,1] .- qc .- qi
    end
    for ilay=1:nlay+1
        for icol=1:ncol
            xc = phys_grid.x_col[icol]
            yc = phys_grid.y_col[icol]
            zc = phys_grid.z[ilay,icol]
            e = phys_grid.el_correspondance[ilay,icol]
            B .= 0.0
            for i=1:mesh.ngl
                for j=1:mesh.ngl
                    for k=1:mesh.ngl
                        ip = mesh.connijk[e,i,j,k]
                        for var=1:6
                            B[var] += A[ip,var]*phys_grid.ω_dyn[ilay,icol,i,j,k]
                        end
                    end
                end
            end
            phys_grid.p[ilay,icol] = B[2]
            phys_grid.ρ[ilay,icol] = B[1]
            phys_grid.t[ilay,icol] = B[3]
            phys_grid.qc[ilay,icol] = B[4] * PhysConst.Mol_mass_water/ PhysConst.Mol_mass_air 
            phys_grid.qi[ilay,icol] = B[5] * PhysConst.Mol_mass_water/ PhysConst.Mol_mass_air
            phys_grid.qv[ilay,icol] = B[6] * PhysConst.Mol_mass_water/ PhysConst.Mol_mass_air
        end
    end
    @info maximum(A[:,2]), minimum(A[:,2]), maximum(phys_grid.p), minimum(phys_grid.p)
    @info maximum(A[:,1]), minimum(A[:,1]), maximum(phys_grid.ρ), minimum(phys_grid.ρ)
    @info maximum(A[:,3]), minimum(A[:,3]), maximum(phys_grid.t), minimum(phys_grid.t)
    @info maximum(A[:,4]), minimum(A[:,4]), maximum(phys_grid.qc), minimum(phys_grid.qc)
    @info maximum(A[:,5]), minimum(A[:,5]), maximum(phys_grid.qi), minimum(phys_grid.qi)
    @info maximum(A[:,6]), minimum(A[:,6]), maximum(phys_grid.qv), minimum(phys_grid.qv)
end

@kernel function interpolate_to_phys_grid_gpu!(x,y,z,connijk,x_c,y_c,z_c,uaux,qe,nlay,ncol,P,Tabs,qc,qi,ρ, P_c,ρ_c,t_c,qc_c,qi_c,qv_c, nvar)
    icol = @index(Group,Linear)
    ilev = @index(Local,Linear)
    e = el_correspondance[ilev,icol]
    T=eltype(x)
    S = zero(T)
    for i=1:ngl
        for j=1:ngl
            for k=1:ngl
                ip = connijk[e,i,j,k]
                S += 1/ (sqrt((x[ip] - x_c)^2 + (y[ip] - y_c)^2 + (z[ip] - z_c)^2))
            end
        end
    end
    for i=1:ngl
        for j=1:ngl
            for k=1:ngl
                ip = connijk[e,i,j,k]
                ω = (1/ (sqrt((x[ip] - x_c)^2 + (y[ip] - y_c)^2 + (z[ip] - z_c)^2)))/S
                P_c[ilev,icol] += P[ip]*ω
                ρ_c[ilev,icol] += ρ*ω
                t_c[ilev,icol] += T*ω
                qc_c[ilev,icol] += qc*ω
                qi_c[ilev,icol] += qi*ω
                qv_c[ilev,icol] += qv*ω
            end
        end
    end


end


function populate_layer_data_from_level_data!(phys_grid,nlay,ncol)

    for ilay =1:nlay
        for icol=1:ncol
            phys_grid.p_lay[ilay,icol] = (phys_grid.p[ilay+1,icol]-phys_grid.p[ilay,icol])/log10(phys_grid.p[ilay+1,icol]/phys_grid.p[ilay,icol]) 
            if (ilay < nlay) 
                phys_grid.t_lay[ilay,icol] = (phys_grid.t[ilay+1,icol] + phys_grid.t[ilay,icol])/2
                phys_grid.qc_lay[ilay,icol] = (phys_grid.qc[ilay+1,icol] + phys_grid.qc[ilay,icol])/2
                phys_grid.qi_lay[ilay,icol] = (phys_grid.qi[ilay+1,icol] + phys_grid.qi[ilay,icol])/2
                phys_grid.qv_lay[ilay,icol] = (phys_grid.qv[ilay+1,icol] + phys_grid.qv[ilay,icol])/2
            else
                phys_grid.t_lay[ilay,icol] = phys_grid.t[ilay+1,icol]
                phys_grid.qc_lay[ilay,icol] = phys_grid.qc[ilay+1,icol]
                phys_grid.qi_lay[ilay,icol] = phys_grid.qi[ilay+1,icol]
                phys_grid.qv_lay[ilay,icol] = phys_grid.qv[ilay+1,icol]
            end
        end
    end
end

@kernel function populate_layer_data_from_level_data!(p,t,qc,qi,qv,p_lay,t_lay,qc_lay,qi_lay,qv_lay,nlay,ncol)
    icol = @index(Group,Linear)
    ilay = @index(Local,Linear)
    
    p[ilay,icol] = (p[ilay+1,icol]-p[ilay,icol])/log10(p[ilay+1,icol]/p[ilay,icol])
    if (ilay < nlay)
        t_lay[ilay,icol] = (t[ilay+1,icol] + t[ilay,icol])/2
        qc_lay[ilay,icol] = (qc[ilay+1,icol] + qc[ilay,icol])/2
        qi_lay[ilay,icol] = (qi[ilay+1,icol] + qi[ilay,icol])/2
        qv_lay[ilay,icol] = (qv[ilay+1,icol] + qv[ilay,icol])/2
    else
        t_lay[ilay,icol] = t[ilay+1,icol]
        qc_lay[ilay,icol] = qc[ilay+1,icol]
        qi_lay[ilay,icol] = qi[ilay+1,icol]
        qv_lay[ilay,icol] = qv[ilay+1,icol]
    end
end

function interpolate_from_phys_grid_cpu!(x,y,z,connijk,phys_grid,flux,flux_interp,nx,ny,ncol,nlay,npoin)

    flux_aux = KernelAbstractions.zeros(CPU(), TFloat, nx,ny,nlay+1)
    for ix=1:nx
        for iy=1:ny
            icol = ix + nx*(iy-1)
            for ilay = 1:nlay+1
                flux_aux[ix,iy,ilay] = flux[ilay,icol]
            end
        end
    end
    
    for ip=1:npoin
        xx=x[ip]
        yy=y[ip]
        zz=z[ip]
        #### find nearest phys_grid column ( this is done to account for warping due to topography the barycentric weights used for interpolation will change from column to column####
        dist = 1000000000.0
        t_col = 1
        for icol=1:ncol
            x_col = phys_grid.x_col[icol]
            y_col = phys_grid.y_col[icol]
            d = sqrt((xx-x_col)^2+(yy-y_col)^2)
            if (d < dist)
                dist = d
                t_col = icol
            end
        end
        Tx = PolynomialInterpolationMatrix(phys_grid.x,phys_grid.ω_x,xx)
        F1 = KernelAbstractions.zeros(CPU(), TFloat, 1, ny,nlay+1)
        for i=1:ny
            for j=1:nlay+1
                F1[:,i,j] .= InterpolateToNewPoints(Tx,flux_aux[:,i,j])
            end
        end
        #@info xx, maximum(Tx), minimum(Tx), maximum(F1), minimum(F1)
        F2 = KernelAbstractions.zeros(CPU(), TFloat, 1, 1, nlay+1)
        Ty = PolynomialInterpolationMatrix(phys_grid.y,phys_grid.ω_y,yy)
        for j=1:nlay+1
            F2[:,:,j] .= InterpolateToNewPoints(Ty,F1[:,:,j])
        end
        Tz = PolynomialInterpolationMatrix(phys_grid.z[:,t_col],phys_grid.ω_z[:,t_col],zz)
        flux_interp[:,ip] .= InterpolateToNewPoints(Tz,F2[:])
    end

end

@kernel function interpolate_from_phys_grid_gpu!(x,y,z,connijk,phys_grid,flux,flux_interp,nx,ny,ncol,nlay)
    ip = @index(Global,Linear)
    T = eltype(x)
    Ti = eltype(ncol)
    xx = x[ip]
    yy = y[ip]
    zz = z[ip]

    dist = T(1000000000.0)
    t_col = Ti(1)

    for icol=1:ncol
        xc = x_c[icol]
        yc = y_c[icol]
        d = sqrt((xx-xc)^2+(yy-yc)^2)
        if (d < dist)
            dist = d
            t_col = icol
        end
    end
    Tx = PolynomialInterpolationMatrix(xp,ω_x,xx)
    for i=1:ny
        for j=1:nlay+1
            F̅ = InterpolateToNewPoints(Tx,flux[:,i,j])
        end
    end
    Ty = PolynomialInterpolationMatrix(yp,ω_y,yy)
    for j=1:nlay+1
        F̃ = InterpolateToNewPoints(Ty,F̅[1,:,j])
    end
    Tz = PolynomialInterpolationMatrix(zp[:,t_col],ω_z[:,t_col],zz)
    flux_interp[ip] = InterpolateToNewPoints(Tz,F̃[1,1,:])
    
end


function construct_atmospheric_state(phys_grid)

end

function compute_radiative_fluxes!(lnew_mesh, mesh, uaux, qe, mp, phys_grid, backend, ::PERT)
    #=if (lnew_mesh)
        if (backend == CPU())
            store_element_maxima!(mesh,mesh.el_max,mesh.el_min,mesh.nelem,mesh.ngl)         
            store_mesh_to_phys_grid_correspondance!(phys_grid,mesh.el_max,mesh.el_min,mesh.nelem,phys_grid.ncol,phys_grid.nlev-1)
            lnew_mesh = false
        else

        end
    end=#
    if (backend == CPU())
        interpolate_to_phys_grid!(mesh,phys_grid,uaux,qe,phys_grid.nlev-1,phys_grid.ncol,@view(uaux[:,end]),mp.Tabs,mp.qc,mp.qi,@view(uaux[:,1]),true)
        
        col_dry = zeros(phys_grid.nlay, phys_grid.ncol)
        rel_hum = zeros(phys_grid.nlay, phys_grid.ncol)
        lon = nothing # This example skips latitude dependent gravity computation
        lat = nothing
        compute_col_gas!(CPU(), phys_grid.p, col_dry, param_set, phys_grid.qv, lat)
        compute_relative_humidity!(CPU(), rel_hum, phys_grid.p_lay, phys_grid.t_lay, param_set, phys_grid.qv)
        
        layerdata = zeros(4, phys_grid.nlay, phys_grid.ncol)
        layerdata[1, :, :] .= col_dry
        layerdata[2, :, :] .= phys_grid.p_lay
        layerdata[3, :, :] .= phys_grid.t_lay
        layerdata[4, :, :] .= rel_hum
        t_sfc = read_data
        as = AtmosphericState(lon, lat, layerdata, phys_grid.p, phys_grid.t, t_sfc, vmr, nothing, nothing),
        sfc_emis,
        sfc_alb,
        cos_zenith,
        irrad,
        bot_at_1,

        flux = zeros(TFloat,phys_grid.nlev,phys_grid.ncol)
        for ilay = 1:phys_grid.nlev
            for icol =1:phys_grid.ncol
                flux[ilay,icol] = icol + ilay
            end
        end

        flux_interp = KernelAbstractions.zeros(backend,TFloat, 1, mesh.npoin)
        interpolate_from_phys_grid_cpu!(mesh.x,mesh.y,mesh.z,mesh.connijk,phys_grid,flux,flux_interp,phys_grid.nx,phys_grid.ny,phys_grid.ncol,phys_grid.nlev-1, mesh.npoin)
        @info maximum(flux), minimum(flux), maximum(flux_interp), minimum(flux_interp)
    else

    end

end

function compute_radiative_fluxes!(lnew_mesh, mesh, uaux, qe, mp, phys_grid, backend, ::TOTAL)
    #=if (lnew_mesh)
        #=if (backend == CPU())
            store_element_maxima!(mesh,mesh.el_max,mesh.el_min,mesh.nelem,mesh.ngl)
            store_mesh_to_phys_grid_correspondance!(phys_grid,mesh.el_max,mesh.el_min,mesh.nelem,phys_grid.ncol,phys_grid.nlev-1)
            lnew_mesh = false
        else

        end=#
    end=#
    if (backend == CPU())
        interpolate_to_phys_grid!(mesh,phys_grid,uaux,qe,phys_grid.nlev-1,phys_grid.ncol,@view(uaux[:,end]),mp.Tabs,mp.qc,mp.qi,@view(uaux[:,1]),false)
        
        flux = rand(TFloat,phys_grid.nlev,phys_grid.ncol)
        
        flux_interp = KernelAbstractions.zeros(backend,TFloat, 1, mesh.npoin)
        interpolate_from_phys_grid_cpu!(mesh.x,mesh.y,mesh.z,mesh.connijk,phys_grid,flux,flux_interp,phys_grid.nx,phys_grid.ny,phys_grid.ncol,phys_grid.nlev-1, mesh.npoin)
    else

    end

end
