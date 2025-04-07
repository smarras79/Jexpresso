include("custom_bcs.jl")

function apply_boundary_conditions!(u, uaux, t,qe,
                                    x, y, z, nx, ny, nz, npoin, npoin_linear, poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ψ, dψ,
                                    xmax, ymax, zmax, xmin, ymin, zmin, RHS, rhs_el, ubdy,
                                    connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_in_elem, bdy_face_type,
                                    connijk, Jef, S_face, S_flux, F_surf, M_surf_inv,
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
                          connijk, Jef, S_face, S_flux, F_surf, M_surf_inv,
                          Tabs, qn,
                          neqs, dirichlet!, neumann, inputs)
    end
end

function apply_periodicity!(u, uaux, t,qe,
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
                           connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_type, RHS, rhs_el,
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
                           connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_type, RHS, rhs_el,
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
                fill!(qbdy, 4325789.0)
                
                user_bc_dirichlet!(@view(uaux[ip,:]), x[ip], y[ip], t, bdy_edge_type[iedge], qbdy, nx_l, ny_l, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])
                
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

    if(inputs[:llaguerre_bc])
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
        end
        for k=1:ngr
            ip = connijk_lag[1,1,k]
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
            ip = connijk_lag[nelem_semi_inf,ngl,k]
            ny_l = 0.0
            nx_l = 1.0
            fill!(qbdy, 4325789.0)     
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
    
    #Map back to u after applying b.c.
    uaux2u!(u, uaux, neqs, npoin)
    
end


function build_custom_bcs_lin_solve!(::NSD_2D, t, x, y, z, nx, ny, nz, npoin, npoin_linear, poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ω,
                                     xmax, ymax, zmax, xmin, ymin, zmin, qbdy, qe,
                                     connijk_lag, bdy_edge_in_elem, bdy_edge_type, RHS, L,
                                     neqs, dirichlet!, neumann, inputs)
    
    for iedge = 1:nedges_bdy

        if (bdy_edge_type[iedge] != "periodic1" &&
            bdy_edge_type[iedge] != "periodic2" &&
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
                           connijk, Jef, S_face, S_flux, F_surf, M_surf_inv,
                           Tabs, qn,
                           neqs, dirichlet!, neumann, inputs)
    #
    # WARNING: Notice that the b.c. are applied to uaux[:,:] and NOT u[:]!
    #          That
    #for ip = 1:npoin
    PhysConst = PhysicalConst{Float64}()
    for iface = 1:nfaces_bdy
        if bdy_face_type[iface] != "periodic1" && bdy_face_type[iface] != "periodic2" && bdy_face_type[iface] != "periodic3" 
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
                #user_bc_neumann!
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
            RHS[:, ieq] .+= S_flux[:,ieq]
        end
    end
    
    #Map back to u after applying b.c.
    uaux2u!(u, uaux, neqs, npoin)
    
end
