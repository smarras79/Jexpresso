include("custom_bcs.jl")

function apply_boundary_conditions!(u, uaux, t,qe,
                                    mesh, metrics, basis,
                                    RHS, rhs_el, ubdy,
                                    ω, SD, neqs, inputs)
  #  build_custom_bcs!(SD, t, mesh, metrics, ω,
  #                    ubdy, uaux, gradu, @view(rhs_el[:,:,:,:]), neqs,
  #                    dirichlet!, neumann,
  #                    zeros(1,1), inputs)

    build_custom_bcs!(SD, t, mesh, metrics, ω,
                      ubdy, uaux, u, qe,
                      @view(RHS[:,:]), @view(rhs_el[:,:,:,:]),
                      neqs, dirichlet!, neumann, inputs)
    
    #end
    
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
                        #@info mesh.x[ip],mesh.y[ip],ieq,t 
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

#=
function yt_build_custom_bcs!(t, mesh, qbdy, q, gradq, bdy_flux, rhs, ::NSD_2D, nvars, metrics, ω, dirichlet!, neumann, L, inputs)
    
    # nedges =  mesh.nedges_int
    # q_size = size(q, 2)
    # L_size = size(L, 1)
    # fill!(bdy_flux, zero(Float64))
        
    for iedge = 1:size(mesh.bdy_edge_in_elem,1)
        iel  = mesh.bdy_edge_in_elem[iedge]

        if mesh.bdy_edge_type[iedge] != "periodic1" && mesh.bdy_edge_type[iedge] != "periodic2"
            tag = mesh.bdy_edge_type[iedge]
            for k=1:mesh.ngl
                ip = mesh.poin_in_bdy_edge[iedge,k]
                ωJacedge = ω[k]*metrics.Jef[iedge,k]
                
                #bdy_flux = zeros(q_size,1)
                x    = mesh.x[ip]
                y    = mesh.y[ip]

                fill!(qbdy, 4325789.0)
                #flags = zeros(q_size,1) #NO, NEVER ALLOCATE INSIDE A LOOP
                if (inputs[:luser_bc]) #in contrast to what?
                    #q[ip,:], flags = dirichlet!(q[ip,:],gradq[:,ip,:],x,y,t,mesh,metrics,tag,qbdy)
                    ipp=1 #ip               
                    qbdy = dirichlet!(@view(q[ip,:]),@view(gradq[:,ipp,:]),x,y,t,mesh,metrics,tag,qbdy,inputs)
                    ##SM change this to set bdy_flux to zero and do not allocate gradq unless neumann is required explicitly by the user
                    bdy_flux .= ωJacedge.*neumann(q[ip,:],gradq[:,ipp,:],x,y,t,mesh,metrics,tag,inputs)
                else
                    q[ip,:] .= 0.0
                end

                
                mm=1; ll=1
                for jj=1:mesh.ngl
                    for ii=1:mesh.ngl
                        if (mesh.connijk[iel,ii,jj] == ip)
                            mm=jj
                            ll=ii
                        end
                    end
                end
                rhs[iel,ll,mm,:] .= rhs[iel,ll,mm,:] .+ bdy_flux[:]
                
                
                for var =1:nvars
                    if !(AlmostEqual(qbdy[var],4325789.0)) # WHAT's this for?
                        @info var,x,y,qbdy[var],
                        
                        rhs[iel,ll,mm,var] = 0.0 #WHAT DOES THIS DO? here is only updated the  `ll` and `mm` row outside of any ll or mm loop
                        
                        q[ip,var]          = qbdy[var]
                    end
                end
                # if (L_size > 1)  ALL THI SABOUT L should exist when we use L
                #, which we don't ever unless we explicitly build it.
                #     for ii=1:mesh.npoin
                #         L[ip,ii] = 0.0
                #     end
                #     L[ip,ip] = 1.0
                # end
                
            end
        end
    end
end
=#
