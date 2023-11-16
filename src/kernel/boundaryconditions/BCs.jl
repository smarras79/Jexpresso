include("custom_bcs.jl")

function apply_boundary_conditions!(u, uaux, t,qe,
                                    mesh, metrics, basis,
                                    RHS, rhs_el, ubdy,
                                    ω, neqs, inputs, SD::NSD_1D)
    if inputs[:lperiodic_1d]
        nothing
    else
        build_custom_bcs!(SD, t, mesh, metrics, ω,
                          ubdy, uaux, u, qe,
                          @view(RHS[:,:]), @view(rhs_el[:,:,:,:]),
                          neqs, dirichlet!, neumann, inputs)
    end
end

function apply_boundary_conditions!(u, uaux, t,qe,
                                    mesh, metrics, basis,
                                    RHS, rhs_el, ubdy,
                                    ω, neqs, inputs, SD::NSD_2D)

    build_custom_bcs!(SD, t, mesh, metrics, ω,
                      ubdy, uaux, u, qe,
                      @view(RHS[:,:]), @view(rhs_el[:,:,:,:]),
                      neqs, dirichlet!, neumann, inputs)
    
end

function apply_periodicity!(u, uaux, t,qe,
                            mesh, metrics, basis,
                            RHS, rhs_el, ubdy,
                            ω, neqs, inputs, SD::NSD_1D)
    nothing
end

function apply_periodicity!(u, uaux, t,qe,
                            mesh, metrics, basis,
                            RHS, rhs_el, ubdy,
                            ω, neqs, inputs, SD::NSD_2D)
    nothing
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

function build_custom_bcs!(::NSD_1D, t, mesh, metrics, ω,
                           qbdy, uaux, u, qe,
                           RHS, rhs_el,
                           neqs, dirichlet!, neumann, inputs)
   ip = 1
   fill!(qbdy, 4325789.0)
   user_bc_dirichlet!(@view(uaux[ip,:]), mesh.x[ip], t, "left", qbdy, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])
   for ieq =1:neqs
      if !AlmostEqual(qbdy[ieq],uaux[ip,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
         uaux[ip,ieq] = qbdy[ieq]
         RHS[ip, ieq] = 0.0
      end
    end
    
    ip=mesh.npoin_linear
    fill!(qbdy, 4325789.0)
    user_bc_dirichlet!(@view(uaux[ip,:]), mesh.x[ip], t, "right", qbdy, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])
    for ieq =1:neqs
      if !AlmostEqual(qbdy[ieq],uaux[ip,ieq]) && !AlmostEqual(qbdy[ieq],4325789.0) # WHAT's this for?
         uaux[ip,ieq] = qbdy[ieq]
         RHS[ip, ieq] = 0.0
      end
    end

    #Map back to u after applying b.c.
    uaux2u!(u, uaux, neqs, mesh.npoin)

end

function build_custom_bcs!(::NSD_2D, t, mesh, metrics, ω,
                           qbdy, uaux, u, qe,
                           RHS, rhs_el,
                           neqs, dirichlet!, neumann, inputs)
    #
    # WARNING: Notice that the b.c. are applied to uaux[:,:] and NOT u[:]!
    #          That
    for iedge = 1:mesh.nedges_bdy 
        iel  = mesh.bdy_edge_in_elem[iedge]
        
        if mesh.bdy_edge_type[iedge] != "periodic1" && mesh.bdy_edge_type[iedge] != "periodic2" && mesh.bdy_edge_type != "Laguerre"
        #if mesh.bdy_edge_type[iedge] == "free_slip"
            
            #tag = mesh.bdy_edge_type[iedge]
            for k=1:mesh.ngl
                ip = mesh.poin_in_bdy_edge[iedge,k]
                nx = metrics.nx[iedge,k]
                ny = metrics.ny[iedge,k]
                fill!(qbdy, 4325789.0)
                #qbdy[:] .= uaux[ip,:]
                #ipp = 1 #ip               
                ###_bc_dirichlet!(qbdy, mesh.x[ip], mesh.y[ip], t, mesh.bdy_edge_type[iedge])

                #dirichlet!(@view(uaux[ip,:]),qbdy, mesh.x[ip], mesh.y[ip], t, metrics.nx[iedge,k], metrics.ny[iedge,k], mesh.bdy_edge_type[iedge], @view(qe[ip,:]), inputs[:SOL_VARS_TYPE]) ###AS IT IS NOW, THIS IS ALLOCATING SHIT TONS. REWRITE to make it with ZERO allocation. hint: It may be due to passing the function but possibly not.
                user_bc_dirichlet!(@view(uaux[ip,:]), mesh.x[ip], mesh.y[ip], t, mesh.bdy_edge_type[iedge], qbdy, nx, ny, @view(qe[ip,:]),inputs[:SOL_VARS_TYPE])
                
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
    
    #Map back to u after applying b.c.
    uaux2u!(u, uaux, neqs, mesh.npoin)
       
end
