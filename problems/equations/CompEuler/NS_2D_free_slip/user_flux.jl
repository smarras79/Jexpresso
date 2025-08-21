function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D,
    q,derivative_y,derivative_x,
    qe::SubArray{Float64},
    mesh::St_mesh,
    ::CL, ::TOTAL; 
    neqs=1, ip=1)   

    F[1] = -q[1]*derivative_y
    G[1] =  q[1]*derivative_x

   #=  ω = params.ω
    ω1 = ω
    ψ = params.basis.ψ
    ψ1 = ψ
    dψ = params.basis.dψ
    dψ1 = dψ
    dξdx = params.metrics.dξdx
    dξdy = params.metrics.dξdy
    dηdx = params.metrics.dηdx
    dηdy = params.metrics.dηdy
    ngl = params.mesh.ngl

    IP = connijk[iel,K,L] =#

    # Exact
   #=  for j = 1:params.mesh.ngl
        for i = 1:params.mesh.ngl# which basis

            J = connijk[iel,i,j]

            dψIJ_dx = dψ[i,K]*ψ1[j,L]*dξdx[iel,K,L] + ψ[i,K]*dψ1[j,L]*dηdx[iel,K,L]
            dψIJ_dy = dψ[i,K]*ψ1[j,L]*dξdy[iel,K,L] + ψ[i,K]*dψ1[j,L]*dηdy[iel,K,L]

            F[1] = F[1] .+ coef[J].*dψIJ_dy.*(-params.uaux[IP])
            G[1] = G[1] .+ coef[J].*dψIJ_dx.*( params.uaux[IP])
                
        end
    end
 =#


    # Inexact
   #=Coef = zeros(TFloat,ngl,ngl)

    for j = 1:ngl
        for i = 1:ngl
            ii = i + (j-1)*ngl
            Coef[i,j] = coef[ii]
        end
    end
      
                dFdξ = 0.0
                dFdη = 0.0

                @turbo for k = 1:ngl
                    dFdξ += dψ[k,K]*Coef[k,L]
                    dFdη += dψ[k,L]*Coef[K,k]
                end

                dξdx_ij = dξdx[iel,K,L]
                dξdy_ij = dξdy[iel,K,L]
                dηdx_ij = dηdx[iel,K,L]
                dηdy_ij = dηdy[iel,K,L]
                
                dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij
               
                dFdy = dFdξ*dξdy_ij + dFdη*dηdy_ij
                
                F[1] = dFdy*(-params.uaux[IP])
                G[1] = dFdx*(params.uaux[IP])
=#
end
