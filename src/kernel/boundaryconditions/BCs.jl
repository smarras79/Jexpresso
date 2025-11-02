include("custom_bcs.jl")

function apply_boundary_conditions_dirichlet!(u, uaux, t,qe,
                                              coords, 
                                              #x, y, z,
                                              nx, ny, nz,
                                              npoin, npoin_linear,
                                              poin_in_bdy_edge, poin_in_bdy_face,
                                              nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ψ, dψ,
                                              xmax, ymax, zmax, xmin, ymin, zmin, RHS, rhs_el, ubdy,
                                              connijk_lag,
                                              bdy_edge_in_elem, bdy_edge_type, bdy_face_in_elem, bdy_face_type,
                                              connijk, Jef, S_face, S_flux, F_surf, M_surf_inv, M_edge_inv, M_inv,
                                              Tabs, qn,
                                              ω, neqs, inputs, AD, SD)
    
    if inputs[:lperiodic_1d] && typeof(SD) == NSD_1D
        
        apply_periodicity!(u, uaux, t,qe,
                           npoin_linear, ψ, dψ,
                           RHS, rhs_el, ubdy,
                           ω, neqs, inputs, AD, SD)
    else
        build_custom_bcs_dirichlet!(SD, t,
                                    coords, 
                                    #x, y, z,
                                    nx, ny, nz, npoin, npoin_linear,
                                    poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ω,
                                    xmax, ymax, zmax, xmin, ymin, zmin, ubdy, uaux, u, qe,
                                    connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_in_elem, bdy_face_type, RHS, rhs_el,
                                    connijk, Jef, S_face, S_flux, F_surf, M_surf_inv, M_edge_inv, M_inv,
                          Tabs, qn,
                          neqs, dirichlet!, neumann, inputs)
    end
end

function apply_boundary_conditions_neumann!(u, uaux, t,qe,
                                            coords, #x, y, z,
                                            nx, ny, nz,
                                            npoin, npoin_linear,
                                            poin_in_bdy_edge, poin_in_bdy_face,
                                            nedges_bdy, nfaces_bdy,
                                            ngl, ngr, nelem_semi_inf, ψ, dψ,
                                            xmax, ymax, zmax, xmin, ymin, zmin,
                                            RHS, rhs_el, ubdy,
                                            connijk_lag, bdy_edge_in_elem,
                                            bdy_edge_type, bdy_face_in_elem, bdy_face_type,
                                            connijk, Jef,
                                            S_face, S_flux, F_surf, M_surf_inv, M_edge_inv, M_inv,
                                            τ_f, wθ,
                                            Tabs, qn,
                                            ω, neqs, inputs, AD, SD)

    build_custom_bcs_neumann!(SD, t,
                              @view(coords[:,:]),
                              #x, y, z,
                              nx, ny, nz, npoin, npoin_linear, 
                              poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ω,
                              xmax, ymax, zmax, xmin, ymin, zmin, ubdy, uaux, u, qe,
                              connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_in_elem, bdy_face_type, RHS, rhs_el,
                              connijk, Jef, S_face, S_flux, F_surf, M_surf_inv, M_edge_inv, M_inv,
                              τ_f, wθ,
                              Tabs, qn,
                              neqs, dirichlet!, neumann, inputs)
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
                                              coords,
                                              nx, ny, nz,
                                              npoin, npoin_linear,
                                              poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy,
                                              ngl, ngr, nelem_semi_inf,
                                              ψ, dψ,
                                              xmax, ymax, zmax, xmin, ymin, zmin,
                                              RHS, rhs_el, ubdy,
                                              connijk_lag, bdy_edge_in_elem, bdy_edge_type,
                                              ω, neqs, inputs, AD, SD)

    if inputs[:lsparse]
        # SM HERE: uncomment this and write it for the Ax=b problem when using Dirichlet.
        build_custom_bcs_lin_solve_sparse!(SD, t, coords, nx, ny, nz, npoin, npoin_linear,
                                           poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy,
                                           ngl, ngr, nelem_semi_inf, ω,
                                           xmax, ymax, zmax, xmin, ymin, zmin, ubdy, qe,
                                           connijk_lag, bdy_edge_in_elem, bdy_edge_type, RHS, L,
                                           neqs, dirichlet!, neumann, inputs)
    else
        # SM HERE: uncomment this and write it for the Ax=b problem when using Dirichlet.
        build_custom_bcs_lin_solve!(SD, t, coords, nx, ny, nz, npoin, npoin_linear,
                                    poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy,
                                    ngl, ngr, nelem_semi_inf, ω,
                                    xmax, ymax, zmax, xmin, ymin, zmin, ubdy, qe,
                                    connijk_lag, bdy_edge_in_elem, bdy_edge_type, RHS, L,
                                    neqs, dirichlet!, neumann, inputs)
    end
    
end

function build_custom_bcs_dirichlet!(::NSD_1D, t,
                                     coords, 
                                     #x, y, z,
                                     nx, ny, nz, npoin, 
                                     npoin_linear, poin_in_bdy_edge, poin_in_bdy_face,
                                     nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ω,
                                     xmax, ymax, zmax, xmin, ymin, zmin, qbdy, uaux, u, qe,
                                     connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_in_elem, bdy_face_type,
                                     RHS, rhs_el,
                                     connijk, Jef, S_face, S_flux, F_surf, M_surf_inv, M_edge_inv, M_inv,
                                     Tabs, qn,
                                     neqs, dirichlet!, neumann, inputs)
    
    ip = 1
    fill!(qbdy, 4325789.0)
    user_bc_dirichlet!(@view(uaux[ip,:]), @view(coords[ip,:]), t, "left", qbdy, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])
    for ieq =1:neqs
        if !AlmostEqual(qbdy[ieq],uaux[ip,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
            uaux[ip,ieq] = qbdy[ieq]
            RHS[ip, ieq] = 0.0
        end
    end
    
    ip=npoin_linear
    fill!(qbdy, 4325789.0)
    user_bc_dirichlet!(@view(uaux[ip,:]), @view(coords[ip,:]), t, "right", qbdy, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])
    for ieq =1:neqs
        if !AlmostEqual(qbdy[ieq],uaux[ip,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
            uaux[ip,ieq] = qbdy[ieq]
            RHS[ip, ieq] = 0.0
        end
    end

    #Map back to u after applying b.c.
    uaux2u!(u, uaux, neqs, npoin)

end

function build_custom_bcs_neumann!(::NSD_1D, t,
                                   coords,
                                   #x, y, z,
                                   nx, ny, nz, npoin, 
                                   npoin_linear, poin_in_bdy_edge, poin_in_bdy_face,
                                   nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ω,
                                   xmax, ymax, zmax, xmin, ymin, zmin, qbdy, uaux, u, qe,
                                   connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_in_elem, bdy_face_type, RHS, rhs_el,
                                   connijk, Jef, S_face, S_flux, F_surf, M_surf_inv, M_edge_inv, M_inv,
                                   τ_f, wθ,
                                   Tabs, qn,
                                   neqs, dirichlet!, neumann, inputs)

    nothing
end

function build_custom_bcs_dirichlet!(::NSD_2D, t,
                                     coords,
                                     #x, y, z,
                                     nx, ny, nz, npoin, 
                                     npoin_linear, poin_in_bdy_edge, poin_in_bdy_face,
                                     nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ω,
                                     xmax, ymax, zmax, xmin, ymin, zmin, qbdy, uaux, u, qe,
                                     connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_in_elem, bdy_face_type, RHS, rhs_el,
                                     connijk, Jef, S_face, S_flux, F_surf, M_surf_inv, M_edge_inv, M_inv,
                                     Tabs, qn,
                                     neqs, dirichlet!, neumann, inputs)

    #
    # WARNING: Notice that the b.c. are applied to uaux[:,:] and NOT u[:]!
    #          That
    for iedge = 1:nedges_bdy 
        iel  = bdy_edge_in_elem[iedge]
        
        if  bdy_edge_type[iedge] != "periodicx" && bdy_edge_type[iedge] != "periodicz" &&
            bdy_edge_type[iedge] != "periodic1" && bdy_edge_type[iedge] != "periodic2" &&
            bdy_edge_type[iedge] != "Laguerre"
            for k=1:ngl
                ip = poin_in_bdy_edge[iedge,k]
                nx_l = nx[iedge,k]
                ny_l = ny[iedge,k]
                fill!(qbdy, 4325789.0)
                
                user_bc_dirichlet!(@view(uaux[ip,:]), @view(coords[ip,:]), t, bdy_edge_type[iedge], qbdy, nx_l, ny_l, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])
                
                for ieq =1:neqs
                    if !AlmostEqual(qbdy[ieq],uaux[ip,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
                        uaux[ip,ieq] = qbdy[ieq]
                        RHS[ip, ieq] = 0.0
                    end
                end
            end
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
                    user_bc_dirichlet!(@view(uaux[ip,:]), @view(coords[ip,:]), t, tag, qbdy, nx_l, ny_l, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])

                    for ieq =1:neqs
                        if !AlmostEqual(qbdy[ieq],uaux[ip,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
                            #@info mesh.x[ip],mesh.y[ip],ieq,qbdy[ieq]
                            uaux[ip,ieq] = qbdy[ieq]
                            RHS[ip, ieq] = 0.0
                        end
                    end
                end
                ip_test = connijk_lag[e,1,1]
                if (coords[ip_test,1] == xmin)
                    for k=1:ngr
                        ip = connijk_lag[e,1,k]
                        ny_l = 0.0
                        nx_l = -1.0
                        fill!(qbdy, 4325789.0)
                        tag = inputs[:laguerre_tag]
                        user_bc_dirichlet!(@view(uaux[ip,:]), @view(coords[ip,:]), t, tag, qbdy, nx_l, ny_l, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])
                        for ieq =1:neqs
                            if !AlmostEqual(qbdy[ieq],uaux[ip,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0)
                                uaux[ip,ieq] = qbdy[ieq]
                                RHS[ip, ieq] = 0.0
                            end
                        end
                    end
                end
                ip_test = connijk_lag[e,ngl,1]
                if (coords[ip_test,1] == xmax)
                    for k =1:ngr
                        ip = connijk_lag[nelem_semi_inf,ngl,k]
                        ny_l = 0.0
                        nx_l = 1.0
                        fill!(qbdy, 4325789.0)
                        tag = inputs[:laguerre_tag]
                        user_bc_dirichlet!(@view(uaux[ip,:]), @view(coords[ip,:]), t, tag, qbdy, nx_l, ny_l, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])

                        for ieq =1:neqs
                            if !AlmostEqual(qbdy[ieq],uaux[ip,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0)
                                uaux[ip,ieq] = qbdy[ieq]
                                RHS[ip, ieq] = 0.0
                            end
                        end
                    end
                end
            end
        end
    end
    uaux2u!(u, uaux, neqs, npoin)
end
function build_custom_bcs_neumann!(::NSD_2D, t,
                                   coords, 
                                   #x, y, z,
                                   nx, ny, nz, npoin,
                                   npoin_linear, poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ω,
                                   xmax, ymax, zmax, xmin, ymin, zmin, qbdy, uaux, u, qe,
                                   connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_in_elem, bdy_face_type, RHS, rhs_el,
                                   connijk, Jef, S_face, S_flux, F_surf, M_surf_inv, M_edge_inv, M_inv,
                                   τ_f, wθ,
                                   Tabs, qn,
                                   neqs, dirichlet!, neumann, inputs)

    # WARNING: Notice that the b.c. are applied to uaux[:,:] and NOT u[:]!
    #          That
    if (inputs[:bdy_fluxes])
        for iedge = 1:nedges_bdy
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
                    bulk_surface_flux!(@view(F_surf[i,:]),
                                       uaux[ip,:], uaux[ip1,:],
                                       qe[ip,:], qe[ip1,:],
                                       θ, θ1,
                                       qn[ip], qn[ip1])
                end

            else
                for i = 1:ngl
                    ip  = poin_in_bdy_edge[iedge,i]
                    e   = bdy_edge_in_elem[iedge]
                    ip1 = connijk[e,i,ngl]
                    
                    user_bc_neumann!(@view(F_surf[i,:]), uaux[ip,:], uaux[ip1,:],
                                     qe[ip,:], qe[ip1,:],
                                     bdy_edge_type[iedge],
                                     @view(coords[ip,:]), inputs[:SOL_VARS_TYPE])
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
end


function build_custom_bcs_lin_solve_sparse!(::NSD_2D, t, coords, nx, ny, nz,
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
                
                user_bc_dirichlet!(@view(RHS[ip,:]),
                                   @view(coords[ip,:]), t,
                                   bdy_edge_type[iedge],
                                   qbdy,
                                   nx_l, ny_l,
                                   @view(qe[ip,:]),
                                   inputs[:SOL_VARS_TYPE])
                
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

    else
        # Replace bdy rows in global L (in sparse storage):
        apply_dirichlet_bc_inplace!(L, poin_in_bdy_edge, ngl)
    end
    
end

function apply_dirichlet_bc_inplace!(L::SparseMatrixCSC, poin_in_bdy_edge, ngl)
 """
    Apply Dirichlet boundary conditions by modifying sparse matrix in-place
    Sets boundary rows to: L[ip, :] = 0, L[ip, ip] = 1
    
    This is efficient when the number of boundary nodes is small relative to matrix size.
    """
    
    # Get boundary node indices
    boundary_nodes = Int[]
    for iedge = 1:size(poin_in_bdy_edge, 1)
        for k = 1:ngl
            ip = poin_in_bdy_edge[iedge, k]
            if ip > 0  # Valid node index
                push!(boundary_nodes, ip)
            end
        end
    end
    
    # Remove duplicates and sort
    boundary_nodes = unique!(sort!(boundary_nodes))
    
    #println("Applying boundary conditions to $(length(boundary_nodes)) nodes")
    
    # Method 1a: Direct modification (works but not most efficient)
    for ip in boundary_nodes
        # Zero out entire row
        for j = 1:size(L, 2)
            L[ip, j] = 0.0
        end
        # Set diagonal to 1
        L[ip, ip] = 1.0
    end
    
    return L
end


function build_custom_bcs_lin_solve!(::NSD_2D, t, coords,
                                     nx, ny, nz,
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
                
                user_bc_dirichlet!(@view(RHS[ip,:]), @view(coords[ip,:]), t,
                                   bdy_edge_type[iedge], qbdy, nx_l, ny_l, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])

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


function build_custom_bcs_dirichlet!(::NSD_3D, t, coords, nx, ny, nz, npoin, npoin_linear, poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ω,
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
                                       @view(coords[ip,:]), 
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
    end

    #Map back to u after applying b.c.
    uaux2u!(u, uaux, neqs, npoin)
end


function build_custom_bcs_neumann!(::NSD_3D, t, coords, nx, ny, nz, npoin, npoin_linear,
                                   poin_in_bdy_edge, poin_in_bdy_face, nedges_bdy, nfaces_bdy, ngl, ngr, nelem_semi_inf, ω,
                                   xmax, ymax, zmax, xmin, ymin, zmin, qbdy, uaux, u, qe,
                                   connijk_lag, bdy_edge_in_elem, bdy_edge_type, bdy_face_in_elem, bdy_face_type, RHS, rhs_el,
                                   connijk, Jef, S_face, S_flux, F_surf, M_surf_inv, M_edge_inv, M_inv,
                                   τ_f, wθ, 
                                   Tabs, qn,
                                   neqs, dirichlet!, neumann, inputs)
    
    PhysConst = PhysicalConst{Float64}() 
    micro = size(Tabs,1)
    for iface = 1:nfaces_bdy
        if (inputs[:bdy_fluxes])
            F_surf .= 0.0
            if (inputs[:bulk_fluxes])
                if (coords[poin_in_bdy_face[iface,3,3],3] == zmin)
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
                #if (coords[poin_in_bdy_face[iface,3,3], 3] == zmin) # FOR YT THIS WOULDN'T WORK WITH TOPOGRAPHY
                    for i = 1:ngl
                        for j = 1:ngl
                            

                            ip  = poin_in_bdy_face[iface,i,j]
                            e   = bdy_face_in_elem[iface]

                            #Inside point
                            iz  = inputs[:ifirst_wall_node_index]
                            ip1 = connijk[e,i,j,iz]
                            
                            θ = 0.0   
                            θ1 = 0.0 
                            qn1 = 0
                            qn2 = 0
                            if (micro > 1)
                                qn1 = qn[ip]
                                qn2 = qn[ip1]
                                if (Tabs[ip] < 1)
                                    θ = uaux[ip,5]/uaux[ip,1]
                                    θ1 = uaux[ip1,5]/uaux[ip1,1]
                                else
                                    θ = Tabs[ip]*(PhysConst.pref/uaux[ip,end])^(1/PhysConst.cpoverR)
                                    θ1 = Tabs[ip1]*(PhysConst.pref/uaux[ip1,end])^(1/PhysConst.cpoverR)
                                end
                            end

                            if (bdy_face_type[iface] == "MOST")
                                
                                ipsfc    = connijk[e,i,j,1]
                                ρ        = uaux[ipsfc, 1]
                                
                                u_inside = uaux[ip1,   2]/ρ
                                v_inside = uaux[ip1,   3]/ρ                                
                                w_inside = uaux[ip1,   4]/ρ

                                vprojectedaux = u_inside*nx[iface,i,j] + v_inside*ny[iface,i,j] + w_inside*nz[iface,i,j]
                                
                                u_inside = u_inside - vprojectedaux*nx[iface,i,j]
                                v_inside = v_inside - vprojectedaux*ny[iface,i,j]
                                w_inside = w_inside - vprojectedaux*nz[iface,i,j]
                                
                                θ_inside = uaux[ip1,   5]/ρ
                                θ_sfc    = uaux[ipsfc, 5]/ρ
                                z_sfc    = coords[ipsfc, 3]
                                                                
                                Δx = coords[ip1, 1] - coords[ipsfc, 1]
                                Δy = coords[ip1, 2] - coords[ipsfc, 2]
                                Δz = coords[ip1, 3] - coords[ipsfc, 3]
                                z_inside = abs(Δx*nx[iface,i,j] + Δy*ny[iface,i,j] + Δz*nz[iface,i,j]) 
                                
                                CM_MOST!(@view(τ_f[iface,i,j,:]), @view(wθ[iface,i,j,1]), ρ, u_inside, v_inside, w_inside, θ_inside, θ_sfc, z_inside)
                                
                                F_surf[i,j,2] = τ_f[iface,i,j,1]
                                F_surf[i,j,3] = τ_f[iface,i,j,2]
                                F_surf[i,j,4] = τ_f[iface,i,j,3]
                                F_surf[i,j,5] = 0.12
                                
                            else
                                user_bc_neumann!(@view(F_surf[i,j,:]), uaux[ip,:], uaux[ip1,:],
                                                 qe[ip,:], qe[ip1,:],
                                                 bdy_face_type[iface],
                                                 @view(coords[ip,:]),
                                                 @view(τ_f[iface,i,j,:]), @view(wθ[iface,i,j,1]), inputs[:SOL_VARS_TYPE], PhysConst;
                                                 θ = θ,
                                                 θ1 = θ1,
                                                 qn0 = qn1,
                                                 qn1=qn2)
                            end
                        end
                    end
                 #end
            end
            compute_surface_integral!(S_face, F_surf, ω, Jef, iface, ngl)            
        end

    end
    #@info maximum(S_face[:,:,:,2]), maximum(S_face[:,:,:,5]), maximum(S_face[:,:,:,6])
    #@info minimum(S_face[:,:,:,2]), minimum(S_face[:,:,:,5]), minimum(S_face[:,:,:,6])
    if (inputs[:bdy_fluxes])
        DSS_surface_integral!(S_flux, S_face, M_surf_inv, nfaces_bdy, ngl,
                              coords[:,3], zmin, connijk, poin_in_bdy_face, bdy_face_in_elem, neqs)
        #@info maximum(S_flux[:,2]), maximum(S_flux[:,5])
        #@info minimum(S_flux[:,2]), minimum(S_flux[:,5])
        for ieq = 1:neqs
            RHS[:, ieq] .+= S_flux[:,ieq] ./ M_inv[:]
        end
    end
end
