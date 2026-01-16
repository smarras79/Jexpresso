include("custom_bcs.jl")

function apply_boundary_conditions!(u, uaux, t,qe,
                                    x, y, z, nx, ny, nz, npoin, npoin_linear, poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ψ, dψ,
                                    xmax, ymax, zmax, xmin, ymin, zmin, RHS, rhs_el, ubdy,
                                    connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_in_elem, bdy_face_type,
                                    connijk, Jef, S_face, S_flux, F_surf, M_surf_inv, M_edge_inv, M_inv,
                                    Tabs, qn,
                                    ω, neqs, inputs, AD, SD)
    
    if inputs[:lperiodic_1d] && typeof(SD) == NSD_1D
       
        apply_periodicity!(u, uaux, t,qe,
                           npoin_linear, ψ, dψ,
                           RHS, rhs_el, ubdy,
                           ω, neqs, inputs, AD, SD)
    else
        build_custom_bcs!(SD, t, x, y, z, nx, ny, nz, npoin, npoin_linear, poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ω,
                          xmax, ymax, zmax, xmin, ymin, zmin, ubdy, uaux, u, qe,
                          connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_in_elem, bdy_face_type, RHS, rhs_el,
                          connijk, Jef, S_face, S_flux, F_surf, M_surf_inv, M_edge_inv, M_inv,
                          Tabs, qn,
                          neqs, dirichlet!, neumann, inputs)
    end
end

function apply_periodicity!(u, uaux, t, qe,
                            npoin_linear, ψ, dψ,
                            RHS, rhs_el, ubdy,
                            ω, neqs, inputs, AD::FD, SD::NSD_1D)

    #this only works for a scalar equation.
    #adjust for systems.
    u[npoin_linear] = u[1]
end


function apply_periodicity!(u, uaux, t,qe,
                            npoin_linear, ψ, dψ,
                            RHS, rhs_el, ubdy,
                            ω, neqs, inputs, AD::ContGal, SD::NSD_1D)
    nothing
end


function apply_periodicity!(u, uaux, t,qe,
                            npoin_linear, ψ, dψ,
                            RHS, rhs_el, ubdy,
                            ω, neqs, inputs, AD::FD, SD::NSD_2D)
    error(" BCs.jl: FD not implemented for NSD_2D yet! ")
end


function apply_periodicity!(u, uaux, t,qe,
                            npoin_linear, ψ, dψ,
                            RHS, rhs_el, ubdy,
                            ω, neqs, inputs, AD::ContGal, SD::NSD_2D)
    nothing
end

function apply_boundary_conditions_lin_solve!(L, t, qe,
                                              x, y, z,
                                              nx, ny, nz,
                                              npoin, npoin_linear,
                                              poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy,
                                              ngl, ngr, nelem_semi_inf,
                                              ψ, dψ,
                                              xmax, ymax, zmax, xmin, ymin, zmin,
                                              RHS, rhs_el, ubdy,
                                              connijk_lag, bdy_edge_in_elem, bdy_edge_type,
                                              ω, neqs, inputs, AD, SD)

    # SM HERE: uncomment this and write it for the Ax=b problem when using Dirichlet.
    build_custom_bcs_lin_solve!(SD, t, x, y, z, nx, ny, nz, npoin, npoin_linear,
                                poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy,
                                ngl, ngr, nelem_semi_inf, ω,
                                xmax, ymax, zmax, xmin, ymin, zmin, ubdy, qe,
                                connijk_lag, bdy_edge_in_elem, bdy_edge_type, RHS, L,
                                neqs, dirichlet!, neumann, inputs)
    
end


function _bc_dirichlet!(qbdy, x, y, t, tag, mesh)

    # WARNING!!!!
    # THIS SHOULD LEVERAGE the bdy node tag rather than checking coordinates
    # REWRITE and make sure that there is no allocation.
    #############
    eps = 10.0
    xmin = mesh.xmin + eps; xmax = mesh.xmax - eps
    ymin = mesh.ymin + eps; ymax = mesh.ymax - eps
    
    #=if ( x <= -4990.0 || x >= 4990.0)
    qbdy[2] = 0.0
    end
    if (y <= 10.0 || y >= 9990.0)
    qbdy[3] = 0.0
    end
    if ((x >= 4990.0 || x <= -4990.0) && (y >= 9990.0 || y <= 10.0))
    qbdy[2] = 0.0
    qbdy[3] = 0.0
    end=#
    if ( x <= xmin || x >= xmax)
        qbdy[2] = 0.0
    end
    if (y <= ymin || y >= ymax)
        qbdy[3] = 0.0
    end
    if ((x >= xmax || x <= xmin) && (y >= ymax || y <= ymin))
        qbdy[2] = 0.0
        qbdy[3] = 0.0
    end
    
end

function build_custom_bcs!(::NSD_1D, t, x, y, z, nx, ny, nz, npoin, npoin_linear, poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ω,
                           xmax, ymax, zmax, xmin, ymin, zmin, qbdy, uaux, u, qe,
                           connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_in_elem, bdy_face_type, RHS, rhs_el,
                           connijk, Jef, S_face, S_flux, F_surf, M_surf_inv, M_edge_inv, M_inv,
                           Tabs, qn,
                           neqs, dirichlet!, neumann, inputs)
    
    ip = 1
    fill!(qbdy, 4325789.0)
    user_bc_dirichlet!(@view(uaux[ip,:]), x[ip], t, "left", qbdy, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])
    for ieq =1:neqs
        if !AlmostEqual(qbdy[ieq],uaux[ip,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
            uaux[ip,ieq] = qbdy[ieq]
            RHS[ip, ieq] = 0.0
        end
    end
    
    ip=npoin_linear
    fill!(qbdy, 4325789.0)
    user_bc_dirichlet!(@view(uaux[ip,:]), x[ip], t, "right", qbdy, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])
    for ieq =1:neqs
        if !AlmostEqual(qbdy[ieq],uaux[ip,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
            uaux[ip,ieq] = qbdy[ieq]
            RHS[ip, ieq] = 0.0
        end
    end

    #Map back to u after applying b.c.
    uaux2u!(u, uaux, neqs, npoin)

end

function build_custom_bcs!(::NSD_2D, t, x, y, z, nx, ny, nz, npoin, npoin_linear, poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ω,
                           xmax, ymax, zmax, xmin, ymin, zmin, qbdy, uaux, u, qe,
                           connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_in_elem, bdy_face_type, RHS, rhs_el,
                           connijk, Jef, S_face, S_flux, F_surf, M_surf_inv, M_edge_inv, M_inv,
                           Tabs, qn,
                           neqs, dirichlet!, neumann, inputs)
    
    # WARNING: Notice that the b.c. are applied to uaux[:,:] and NOT u[:]!
    #          That
    for iedge = 1:nedges_bdy 
        iel  = bdy_edge_in_elem[iedge]
        
        if  bdy_edge_type[iedge] != "periodicx" && bdy_edge_type[iedge] != "periodicz" &&
            bdy_edge_type[iedge] != "periodic1" && bdy_edge_type[iedge] != "periodic2" &&
            bdy_edge_type[iedge] != "Laguerre"
            #if mesh.bdy_edge_type[iedge] == "free_slip"
            
            #tag = mesh.bdy_edge_type[iedge]
            for k=1:ngl
                ip = poin_in_bdy_edge[iedge,k]
                nx_l = nx[iedge,k]
                ny_l = ny[iedge,k]
                fill!(qbdy, 4325789.0)
                
                user_bc_dirichlet!(@view(uaux[ip,:]), x[ip], y[ip], t, bdy_edge_type[iedge], qbdy, nx_l, ny_l, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])
                
                for ieq =1:neqs
                    if !AlmostEqual(qbdy[ieq],uaux[ip,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
                        #@info x[ip],y[ip],ieq,qbdy[ieq] 
                        uaux[ip,ieq] = qbdy[ieq]
                        RHS[ip, ieq] = 0.0
                    end
                end
            end
        end
        if (inputs[:bdy_fluxes])
            F_surf .= 0.0
            if (inputs[:bulk_fluxes])
                for i = 1:ngl
                    ip  = poin_in_bdy_edge[iedge,i]
                    e   = bdy_edge_in_elem[iedge]
                    ip1 = connijk[e,i,2]
                    if (Tabs[ip] < 1)
                        θ = 0.0
                        θ1 = 0.0
                    else
                        θ = Tabs[ip]*(PhysConst.pref/uaux[ip,end])^(1/PhysConst.cpoverR)
                        θ1 = Tabs[ip1]*(PhysConst.pref/uaux[ip1,end])^(1/PhysConst.cpoverR)
                    end
                    bulk_surface_flux!(@view(F_surf[i,:]), uaux[ip,:], uaux[ip1,:], qe[ip,:], qe[ip1,:], θ, θ1, qn[ip], qn[ip1])
                end

            else
                for i = 1:ngl
                    ip  = poin_in_bdy_edge[iedge,i]
                    e   = bdy_edge_in_elem[iedge]
                    ip1 = connijk[e,i,2]

                    user_bc_neumann!(@view(F_surf[i,:]), uaux[ip,:], uaux[ip1,:], qe[ip,:], qe[ip1,:], bdy_edge_type[iedge], x[ip], y[ip], inputs[:SOL_VARS_TYPE])
                end
            end
            compute_segment_integral!(S_face, F_surf, ω, Jef, iedge, ngl)
        end
    end
    if (inputs[:bdy_fluxes])
        DSS_segment_integral!(S_flux, S_face, M_edge_inv, nedges_bdy, ngl, connijk, poin_in_bdy_edge, bdy_edge_in_elem)
        #@info maximum(S_flux[:,2]), maximum(S_flux[:,5]), maximum(S_flux[:,6])
        #@info minimum(S_flux[:,2]), minimum(S_flux[:,5]), minimum(S_flux[:,6])
        for ieq = 1:neqs
            RHS[:, ieq] .+= S_flux[:,ieq] ./ M_inv[:]
        end
    
    end

    if(inputs[:llaguerre_bc])
        if (nelem_semi_inf >0)
            for e=1:nelem_semi_inf
                for i=1:ngl
                    ip = connijk_lag[e,i,ngr]
                    ny_l = 1.0
                    nx_l = 0.0
                    fill!(qbdy, 4325789.0)
                    tag = inputs[:laguerre_tag]
                    user_bc_dirichlet!(@view(uaux[ip,:]), x[ip], y[ip], t, tag, qbdy, nx_l, ny_l, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])
                
                    for ieq =1:neqs
                        if !AlmostEqual(qbdy[ieq],uaux[ip,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
                            #@info mesh.x[ip],mesh.y[ip],ieq,qbdy[ieq]
                            uaux[ip,ieq] = qbdy[ieq]
                            RHS[ip, ieq] = 0.0
                        end
                    end
                end
                ip_test = connijk_lag[e,1,1]
                if (x[ip_test] == xmin)
                    for k=1:ngr
                        ip = connijk_lag[e,1,k]
                        ny_l = 0.0
                        nx_l = -1.0
                        fill!(qbdy, 4325789.0)
                        tag = inputs[:laguerre_tag]
                        user_bc_dirichlet!(@view(uaux[ip,:]), x[ip], y[ip], t, tag, qbdy, nx_l, ny_l, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])
            
                        for ieq =1:neqs
                            if !AlmostEqual(qbdy[ieq],uaux[ip,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
                    #@info mesh.x[ip],mesh.y[ip],ieq,qbdy[ieq]
                                uaux[ip,ieq] = qbdy[ieq]
                                RHS[ip, ieq] = 0.0
                            end
                        end
                    end
                end
                ip_test = connijk_lag[e,ngl,1]
                if (x[ip_test] == xmax)
                    for k =1:ngr 
                        ip = connijk_lag[nelem_semi_inf,ngl,k]
                        ny_l = 0.0
                        nx_l = 1.0
                        fill!(qbdy, 4325789.0)     
                        tag = inputs[:laguerre_tag]
                        user_bc_dirichlet!(@view(uaux[ip,:]), x[ip], y[ip], t, tag, qbdy, nx_l, ny_l, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])
            
                        for ieq =1:neqs
                            if !AlmostEqual(qbdy[ieq],uaux[ip,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
                    #@info mesh.x[ip],mesh.y[ip],ieq,qbdy[ieq]
                                uaux[ip,ieq] = qbdy[ieq]
                                RHS[ip, ieq] = 0.0
                            end
                        end
                    end
                end
            end
        end
    end
    
    #Map back to u after applying b.c.
    uaux2u!(u, uaux, neqs, npoin)
    
end


function build_custom_bcs_lin_solve!(::NSD_2D, t, x, y, z, nx, ny, nz,
                                     npoin, npoin_linear, poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy,
                                     ngl, ngr, nelem_semi_inf, ω,
                                     xmax, ymax, zmax, xmin, ymin, zmin, qbdy, qe,
                                     connijk_lag, bdy_edge_in_elem, bdy_edge_type, RHS, L,
                                     neqs, dirichlet!, neumann, inputs)
    
    for iedge = 1:nedges_bdy

        if (bdy_edge_type[iedge] != "periodicx" && bdy_edge_type[iedge] != "periodic1" &&
            bdy_edge_type[iedge] != "periodicz" && bdy_edge_type[iedge] != "periodic3" &&
            bdy_edge_type[iedge] != "Laguerre")
            for k=1:ngl
                ip = poin_in_bdy_edge[iedge,k]
                nx_l = nx[iedge,k]
                ny_l = ny[iedge,k]
                fill!(qbdy, 4325789.0)
                
                user_bc_dirichlet!(@view(RHS[ip,:]), x[ip], y[ip], t, bdy_edge_type[iedge], qbdy, nx_l, ny_l, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])

                for ip1 = 1:npoin
                    L[ip,ip1] = 0.0
                end
                L[ip,ip] = 1.0
                for ieq=1:neqs
                    RHS[ip,ieq] = qbdy[ieq]
                end

            end
        end
    end
#@mystop(BCs.jl)
    if ("Laguerre" in bdy_edge_type)
        for k=1:ngr
            ip = connijk_lag[1, 1, k]
            for ip1 = 1:npoin
                L[ip,ip1] = 0.0
            end
            L[ip,ip] = 1.0
            RHS[ip] = 0.0
        end

        for k=1:ngr
            ip = connijk_lag[nelem_semi_inf, ngl, k]
            for ip1 = 1:npoin
                L[ip,ip1] = 0.0
            end
            L[ip,ip] = 1.0
            RHS[ip] = 0.0
        end
        
        for e=1:nelem_semi_inf
            for i=1:ngl
                ip = connijk_lag[e, i, ngr]
                for ip1 = 1:npoin
                    L[ip,ip1] = 0.0
                end
                L[ip,ip] = 1.0
                RHS[ip] = 0.0
            end
        end
    end
    
end


function build_custom_bcs!(::NSD_3D, t, x, y, z, nx, ny, nz, npoin, npoin_linear, poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ω,
                           xmax, ymax, zmax, xmin, ymin, zmin, qbdy, uaux, u, qe,
                           connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_in_elem, bdy_face_type, RHS, rhs_el,
                           connijk, Jef, S_face, S_flux, F_surf, M_surf_inv, M_edge_inv, M_inv,
                           Tabs, qn,
                           neqs, dirichlet!, neumann, inputs)
    #
    # WARNING: Notice that the b.c. are applied to uaux[:,:] and NOT u[:]!
    #          That
    #for ip = 1:npoin
    PhysConst = PhysicalConst{Float64}()
    for iface = 1:nfaces_bdy
        if (bdy_face_type[iface] != "periodicx" && bdy_face_type[iface] != "periodic1" &&
            bdy_face_type[iface] != "periodicz" && bdy_face_type[iface] != "periodic2" &&
            bdy_face_type[iface] != "periodicy" && bdy_face_type[iface] != "periodic3" )
            for i=1:ngl
                for j=1:ngl
                    fill!(qbdy, 4325789.0)
                    ip = poin_in_bdy_face[iface,i,j]
                    user_bc_dirichlet!(@view(uaux[ip,:]),
                                       x[ip], y[ip], z[ip],
                                       t, bdy_face_type[iface], qbdy,
                                       nx[iface,i,j], ny[iface,i,j], nz[iface,i,j],
                                       xmin, xmax,
                                       ymin, ymax,
                                       zmin, zmax,
                                       @view(qe[ip,:]), inputs[:SOL_VARS_TYPE])
                    for ieq =1:neqs
                        if !AlmostEqual(qbdy[ieq],uaux[ip,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
                            uaux[ip,ieq] = qbdy[ieq]
                            RHS[ip, ieq] = 0.0
                        end
                    end
                end
            end
        end
        if (inputs[:bdy_fluxes])
            F_surf .= 0.0
            if (inputs[:bulk_fluxes])
                if (z[poin_in_bdy_face[iface,3,3]] == zmin)
                    for i = 1:ngl
                        for j = 1:ngl
                            ip  = poin_in_bdy_face[iface,i,j]
                            e   = bdy_face_in_elem[iface]
                            ip1 = connijk[e,i,j,2]
                            if (Tabs[ip] < 1)
                                θ = 0.0
                                θ1 = 0.0
                            else
                                θ = Tabs[ip]*(PhysConst.pref/uaux[ip,end])^(1/PhysConst.cpoverR)
                                θ1 = Tabs[ip1]*(PhysConst.pref/uaux[ip1,end])^(1/PhysConst.cpoverR)
                            end
                            bulk_surface_flux!(@view(F_surf[i,j,:]), uaux[ip,:], uaux[ip1,:], qe[ip,:], qe[ip1,:], θ, θ1, qn[ip], qn[ip1])                        
                            #@info F_surf[i,j,:]
                        end
                    end
                end
            
            else
                if (z[poin_in_bdy_face[iface,3,3]] == zmin)
                    for i = 1:ngl
                        for j = 1:ngl
                            ip  = poin_in_bdy_face[iface,i,j]
                            e   = bdy_face_in_elem[iface]
                            ip1 = connijk[e,i,j,2]
                            user_bc_neumann!(@view(F_surf[i,j,:]), uaux[ip,:], uaux[ip1,:], qe[ip,:], qe[ip1,:], bdy_face_type[iface], x[ip], y[ip], z[ip], inputs[:SOL_VARS_TYPE])
                        end
                    end
                 end
            end
            #@info F_surf
            compute_surface_integral!(S_face, F_surf, ω, Jef, iface, ngl)            
        end

    end
    #@info maximum(S_face[:,:,:,2]), maximum(S_face[:,:,:,5]), maximum(S_face[:,:,:,6])
    #@info minimum(S_face[:,:,:,2]), minimum(S_face[:,:,:,5]), minimum(S_face[:,:,:,6])
    if (inputs[:bdy_fluxes])
        DSS_surface_integral!(S_flux, S_face, M_surf_inv, nfaces_bdy, ngl, z, zmin, connijk, poin_in_bdy_face, bdy_face_in_elem)
        #@info maximum(S_flux[:,2]), maximum(S_flux[:,5]), maximum(S_flux[:,6])
        #@info minimum(S_flux[:,2]), minimum(S_flux[:,5]), minimum(S_flux[:,6])
        for ieq = 1:neqs
            RHS[:, ieq] .+= S_flux[:,ieq] ./ M_inv[:]
        end
    end
    
    #Map back to u after applying b.c.
    uaux2u!(u, uaux, neqs, npoin)
    
end


#  Vlasov-Poisson B.C. (no laguerre yet)

# 1D:
# spatial
function build_custom_bcs_VP!(::NSD_1D, t, x, y, z, nx, ny, nz, npoin, npoin_extra, npoin_linear, npoin_linear_extra, poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy, nelem, nelem_extra, ngl, ngl_extra, ngr, nelem_semi_inf, ω,
                           xmax, ymax, zmax, xmin, ymin, zmin, qbdy, uaux, u, qe,
                           connijk, connijk_extra, connijk_lag, 
                           bdy_edge_in_elem, bdy_edge_type, bdy_face_type, RHS, rhs_el,
                           neqs, dirichlet!, neumann, inputs)
    ip = 1
    fill!(qbdy, 4325789.0)
    for iel_extra = 1:nelem_extra
        for iextra = 1:npoin_extra
            ip_extra = connijk_extra[iel_extra, iextra, 1]
            ind1 = (ip-1)*npoin_extra + ip_extra
            user_bc_dirichlet!(@view(uaux[ind1,:]), t, "left", qbdy, @view(qe[ind1,:]),inputs)
            for ieq =1:neqs
                if !AlmostEqual(qbdy[ieq],uaux[ind1,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
                    uaux[ind1,ieq] = qbdy[ieq]
                    RHS[ind1, ieq] = 0.0
                end
            end
        end
    end
    
    ip = npoin_linear
    fill!(qbdy, 4325789.0)
    for iel_extra = 1:nelem_extra
        for iextra = 1:npoin_extra
            ip_extra = connijk_extra[iel_extra, iextra, 1]
            ind2 = (ip-1)*npoin_extra + ip_extra
            user_bc_dirichlet!(@view(uaux[ind2,:]), t, "right", qbdy, @view(qe[ind2,:]),inputs)
            for ieq =1:neqs
                if !AlmostEqual(qbdy[ieq],uaux[ind2,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
                    uaux[ind2,ieq] = qbdy[ieq]
                    RHS[ind2, ieq] = 0.0
                end
            end
        end
    end

    #Map back to u after applying b.c.
    uaux2u!(u, uaux, neqs, npoin, npoin_extra)

end
# velocity
function build_custom_bcs_VP_extra!(::NSD_1D, t, x, y, z, nx, ny, nz, npoin, npoin_extra, npoin_linear, npoin_linear_extra, poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy, nelem, nelem_extra, ngl, ngl_extra, ngr, nelem_semi_inf, ω,
                           xmax, ymax, zmax, xmin, ymin, zmin, qbdy, uaux, u, qe,
                           connijk, connijk_extra, connijk_lag, 
                           bdy_edge_in_elem, bdy_edge_type, bdy_face_type, RHS, rhs_el,
                           neqs, dirichlet!, neumann, inputs)
    ip_extra = 1
    fill!(qbdy, 4325789.0)
    for iel = 1:nelem
        for i = 1:ngl
            ip = connijk[iel,i,1]
            ind1 = (ip-1)*npoin_extra + ip_extra
            user_bc_dirichlet_extra!(@view(uaux[ind1,:]), t, "left", qbdy, @view(qe[ind1,:]),inputs)
            for ieq =1:neqs
                if !AlmostEqual(qbdy[ieq],uaux[ind1,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
                    uaux[ind1,ieq] = qbdy[ieq]
                    RHS[ind1, ieq] = 0.0
                end
            end
        end
    end
    
    ip_extra = npoin_linear_extra
    fill!(qbdy, 4325789.0)
    for iel = 1:nelem
        for i = 1:ngl
            ip = connijk[iel,i,1]
            ind2 = (ip-1)*npoin_extra + ip_extra
            user_bc_dirichlet_extra!(@view(uaux[ind2,:]), t, "right", qbdy, @view(qe[ind2,:]),inputs)
            for ieq =1:neqs
                if !AlmostEqual(qbdy[ieq],uaux[ind2,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
                    uaux[ind2,ieq] = qbdy[ieq]
                    RHS[ind2, ieq] = 0.0
                end
            end
        end
    end

    #Map back to u after applying b.c.
    uaux2u!(u, uaux, neqs, npoin, npoin_extra)

end

# 2D:
# spatial
function build_custom_bcs_VP!(::NSD_2D, t, x, y, z, nx, ny, nz,
                           vx, vy, vz, vnx, vny, vnz, 
                           npoin, npoin_extra, npoin_linear, npoin_linear_extra, 
                           poin_in_bdy_edge, poin_in_bdy_edge_extra, poin_in_bdy_face, poin_in_bdy_face_extra, 
                           nedges_bdy, nedges_bdy_extra, nfaces_bdy, ngl, ngl_extra, ngr, ngr_extra, nelem_extra, nelem_semi_inf, ω, ω_extra,
                           xmax, ymax, zmax, xmin, ymin, zmin, qbdy, uaux, u, qe,
                           connijk, connijk_extra, connijk_lag, 
                           bdy_edge_in_elem, bdy_edge_type, bdy_face_type,
                           bdy_edge_in_elem_extra, bdy_edge_type_extra, bdy_face_type_extra, 
                           RHS, rhs_el,
                           neqs, dirichlet!, neumann, inputs)
    #
    # WARNING: Notice that the b.c. are applied to uaux[:,:] and NOT u[:]!
    #          That
    for iedge = 1:nedges_bdy 
        iel  = bdy_edge_in_elem[iedge]
        
        if bdy_edge_type[iedge] != "periodic1" && bdy_edge_type[iedge] != "periodic2" && bdy_edge_type[iedge] != "Laguerre"
            #if mesh.bdy_edge_type[iedge] == "free_slip"
            
            #tag = mesh.bdy_edge_type[iedge]
            for k=1:ngl
                ip = poin_in_bdy_edge[iedge,k]
                nx_l = nx[iedge,k]
                ny_l = ny[iedge,k]

                for iel_extra = 1:nelem_extra
                    for i = 1:ngl_extra
                        for j = 1:ngl_extra
                            ip_extra = connijk_extra[iel_extra,i,j]

                            ind = (ip-1)*npoin_extra + ip_extra
                            fill!(qbdy, 4325789.0)
                            
                            user_bc_dirichlet!(@view(uaux[ind,:]), x[ip], y[ip], t, bdy_edge_type[iedge], qbdy, nx_l, ny_l, @view(qe[ind,:]),inputs[:SOL_VARS_TYPE])
                            
                            for ieq =1:neqs
                                if !AlmostEqual(qbdy[ieq],uaux[ind,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
                                    #@info mesh.x[ip],mesh.y[ip],ieq,qbdy[ieq] 
                                    uaux[ind,ieq] = qbdy[ieq]
                                    RHS[ind, ieq] = 0.0
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    #Map back to u after applying b.c.
    uaux2u!(u, uaux, neqs, npoin, npoin_extra)
    
end
# velocity
function build_custom_bcs_VP_extra!(::NSD_2D, t, x, y, z, nx, ny, nz,
                           vx, vy, vz, vnx, vny, vnz, 
                           npoin, npoin_extra, npoin_linear, npoin_linear_extra, 
                           poin_in_bdy_edge, poin_in_bdy_edge_extra, poin_in_bdy_face, poin_in_bdy_face_extra, 
                           nedges_bdy, nedges_bdy_extra, nfaces_bdy, ngl, ngl_extra, ngr, ngr_extra, nelem, nelem_semi_inf, ω, ω_extra,
                           xmax, ymax, zmax, xmin, ymin, zmin, qbdy, uaux, u, qe,
                           connijk, connijk_extra, connijk_lag, 
                           bdy_edge_in_elem, bdy_edge_type, bdy_face_type,
                           bdy_edge_in_elem_extra, bdy_edge_type_extra, bdy_face_type_extra, 
                           RHS, rhs_el,
                           neqs, dirichlet!, neumann, inputs)

    for iedge_extra = 1:nedges_bdy_extra 
        iel_extra  = bdy_edge_in_elem_extra[iedge]
        
        if bdy_edge_type_extra[iedge_extra] != "periodic1" && bdy_edge_type_extra[iedge_extra] != "periodic2" && bdy_edge_type_extra[iedge_extra] != "Laguerre"
            
            for k=1:ngl_extra
                ip_extra = poin_in_bdy_edge_extra[iedge_extra,k]
                nx_l_extra = nx_extra[iedge_extra,k]
                ny_l_extra = ny_extra[iedge_extra,k]

                for iel = 1:nelem
                    for i = 1:ngl
                        for j = 1:ngl
                            ip = connijk[iel,i,j]

                            ind = (ip-1)*npoin_extra + ip_extra
                            fill!(qbdy, 4325789.0)
                            
                            user_bc_dirichlet_extra!(@view(uaux[ind,:]), vx[ip], vy[ip], t, bdy_edge_type_extra[iedge_extra], qbdy, nx_l_extra, ny_l_extra, @view(qe[ind,:]),inputs[:SOL_VARS_TYPE])
                            
                            for ieq =1:neqs
                                if !AlmostEqual(qbdy[ieq],uaux[ind,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0)
                                    uaux[ind,ieq] = qbdy[ieq]
                                    RHS[ind, ieq] = 0.0
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    #Map back to u after applying b.c.
    uaux2u!(u, uaux, neqs, npoin, npoin_extra)
    
end

# 3D:
# spatial
function build_custom_bcs_VP!(::NSD_3D, t, x, y, z, nx, ny, nz, 
                             npoin, npoin_extra, npoin_linear, npoin_linear_extra,
                             poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy, 
                             ngl, ngl_extra, ngr, ngr_extra, nelem_extra, nelem_semi_inf, ω,
                             xmax, ymax, zmax, xmin, ymin, zmin, qbdy, uaux, u, qe,
                             connijk, connijk_extra, connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_type, RHS, rhs_el,
                             neqs, dirichlet!, neumann, inputs)
    #
    # WARNING: Notice that the b.c. are applied to uaux[:,:] and NOT u[:]!
    for iface = 1:nfaces_bdy
        if (bdy_face_type[iface] != "periodicx" && bdy_face_type[iface] != "periodic1" &&
            bdy_face_type[iface] != "periodicz" && bdy_face_type[iface] != "periodic2" &&
            bdy_face_type[iface] != "periodicy" && bdy_face_type[iface] != "periodic3" )
            for i=1:ngl
                for j=1:ngl
                    fill!(qbdy, 4325789.0)
                    ip = poin_in_bdy_face[iface,i,j]

                    for iel_extra = 1:nelem_extra
                        for n = 1:ngl_extra,m = 1:ngl_extra,l = 1:ngl_extra
                            ip_extra = connijk_extra[iel_extra,l,m,n]

                            ind = (ip-1)*npoin_extra + ip_extra

                            user_bc_dirichlet!(@view(uaux[ind,:]),
                                            x[ip], y[ip], z[ip],
                                            t, bdy_face_type[iface],  qbdy,
                                            nx[iface,i,j], ny[iface,i,j], nz[iface,i,j],
                                            xmin, xmax,
                                            ymin, ymax,
                                            zmin, zmax,
                                            @view(qe[ind,:]), inputs[:SOL_VARS_TYPE])
                        
                            for ieq =1:neqs
                                if !AlmostEqual(qbdy[ieq],uaux[ind,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
                                    uaux[ind,ieq] = qbdy[ieq]
                                    RHS[ind, ieq] = 0.0
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    #Map back to u after applying b.c.
    uaux2u!(u, uaux, neqs, npoin, npoin_extra)
    
end
# velocity
function build_custom_bcs_VP_extra!(::NSD_3D, t, x, y, z, nx, ny, nz, 
                             vx, vy, vz, vnx, vny, vnz, 
                             npoin, npoin_extra, npoin_linear, npoin_linear_extra,
                             poin_in_bdy_edge, poin_in_bdy_face_extra, nedges_bdy, nfaces_bdy_extra, 
                             ngl, ngl_extra, ngr, ngr_extra, nelem, nelem_semi_inf, ω_extra,
                             vxmax, vymax, vzmax, vxmin, vymin, vzmin, qbdy, uaux, u, qe,
                             connijk, connijk_extra, connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_type_extra, RHS, rhs_el,
                             neqs, dirichlet!, neumann, inputs)
    #
    # WARNING: Notice that the b.c. are applied to uaux[:,:] and NOT u[:]!
    for iface = 1:nfaces_bdy_extra
        if (bdy_face_type_extra[iface] != "periodicx" && bdy_face_type_extra[iface] != "periodic1" &&
            bdy_face_type_extra[iface] != "periodicz" && bdy_face_type_extra[iface] != "periodic2" &&
            bdy_face_type_extra[iface] != "periodicy" && bdy_face_type_extra[iface] != "periodic3" )
            for i=1:ngl_extra
                for j=1:ngl_extra
                    fill!(qbdy, 4325789.0)
                    ip_extra = poin_in_bdy_face_extra[iface,i,j]

                    for iel = 1:nelem
                        for n = 1:ngl,m = 1:ngl,l = 1:ngl
                            ip = connijk[iel,l,m,n]

                            ind = (ip-1)*npoin_extra + ip_extra

                            user_bc_dirichlet_extra!(@view(uaux[ind,:]),
                                            vx[ip_extra], vy[ip_extra], vz[ip_extra],
                                            t, bdy_face_type_extra[iface],  qbdy,
                                            vnx[iface,i,j], vny[iface,i,j], vnz[iface,i,j],
                                            vxmin, vxmax,
                                            vymin, vymax,
                                            vzmin, vzmax,
                                            @view(qe[ind,:]), inputs[:SOL_VARS_TYPE])
                        
                            for ieq =1:neqs
                                if !AlmostEqual(qbdy[ieq],uaux[ind,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
                                    uaux[ind,ieq] = qbdy[ieq]
                                    RHS[ind, ieq] = 0.0
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    #Map back to u after applying b.c.
    uaux2u!(u, uaux, neqs, npoin, npoin_extra)
    
end