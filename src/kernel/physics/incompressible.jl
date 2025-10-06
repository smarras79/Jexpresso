using SparseArrays

function draw_results(params, Coef, time, result)

    SD      = params.SD
    neqs    = params.neqs
    ngl     = params.mesh.ngl
    nelem   = params.mesh.nelem
    ω       = params.ω
    ω1      = ω
    ψ       = params.basis.ψ
    ψ1      = ψ
    dψ      = params.basis.dψ
    dψ1     = dψ
    dξdx    = params.metrics.dξdx
    dξdy    = params.metrics.dξdy
    dηdx    = params.metrics.dηdx
    dηdy    = params.metrics.dηdy
    connijk = params.mesh.connijk 
    x = params.mesh.x
    y = params.mesh.y
    xmin = params.xmin
    xmax = params.xmax
    ymin = params.ymin
    ymax = params.ymax

    inputs = params.inputs

    velocity = zeros(TFloat,maximum(connijk))
    velocity2 = zeros(TFloat,maximum(connijk))
    count = zeros(TFloat,maximum(connijk))
    for iel = 1:nelem
        for j = 1:ngl, i = 1:ngl # points

            ip = connijk[iel,i,j]

            for n = 1:ngl, m = 1:ngl # basis functions

                J = connijk[iel,m,n]
                dψij_dy = dψ[m,i]*ψ1[n,j]*dξdy[iel,i,j] + ψ[m,i]*dψ1[n,j]*dηdy[iel,i,j]
                dψij_dx = dψ[m,i]*ψ1[n,j]*dξdx[iel,i,j] + ψ[m,i]*dψ1[n,j]*dηdx[iel,i,j]
                velocity[ip] = velocity[ip] - Coef[J]*dψij_dy
                velocity2[ip] = velocity2[ip] + Coef[J]*dψij_dx

            end
            count[ip] = count[ip] + 1
        end
    end
    velocity = velocity./count
    velocity2 = velocity2./count

    outvarsref = Array{Union{Nothing, String}}(nothing, neqs)
    qvars = ("ω")
    for i = 1:length(outvarsref)
        outvarsref[i] = string(qvars[i], "_ref")
    end
    
    print(time)
    print(" | ")
    print(maximum(velocity))
    print(" | ")
    println(maximum(Coef))

    if (result == "velocity")
        data = velocity
    else
        data = Coef
    end

    N = 2
    res = zeros(TFloat, 3, maximum(connijk))

    if (round(time,digits=3) == inputs[:tend])

        for iel = 1:nelem
            for j = 1:ngl, i = 1:ngl # points
                ip = connijk[iel,i,j]
                if (round(x[ip],digits=4) == (xmin+xmax)/2)
                    res[1,ip] = y[ip]
                    res[2,ip] = velocity[ip]
                    res[3,ip] = velocity2[ip]
                end
            end
        end

        ind = sortperm(reshape(res[1,:],(1,maximum(connijk))), dims = 2)
        println("================================================")
        println(res[2,ind])
        println(res[3,ind])
        println("================================================")

        write_vtk_ref(SD, params.mesh, data, string(time), inputs[:output_dir]; nvar=length(params.qp.qe[1,:]), outvarsref=outvarsref)
    end

end

function resetRHS_pToZero_viscous!(params, SD::NSD_2D)
    fill!(@view(params.RHS_laplacian_p[:,:,:]),  zero(params.T))
    fill!(@view(params.rhs_laplacian_el_p[:,:,:,:]), zero(params.T))
    fill!(@view(params.rhs_el_p[:,:,:]), zero(params.T))
    fill!(@view(params.RHS_p[:,:]),     zero(params.T))
end

function resetSeg_pToZero_viscous!(params, SD::NSD_2D)
    fill!(@view(params.seg_p[:]),  zero(params.T))
end

function compute_∇ψ!(params, connijk, x, y, F_data)

    ngl = params.mesh.ngl
    nelem = params.mesh.nelem
    npoin = params.mesh.npoin

    ψ = params.basis.ψ
    dψ = params.basis.dψ
    dξdx = params.metrics.dξdx
    dξdy = params.metrics.dξdy
    dηdx = params.metrics.dηdx
    dηdy = params.metrics.dηdy

    poisson(params, connijk, x, y, 
            @view(params.poisson[:,:]), 
            @view(params.RHS_laplacian_p[:,:,:]), @view(params.RHS_p[:,:]))

    # compute derivative of stream function
    fill!(F_data,zero(params.T))
    count = zeros(TFloat, npoin)

    for iel = 1:nelem
        for l = 1:ngl,k = 1:ngl
            ip = connijk[iel,k,l]

            for j = 1:params.mesh.ngl
                for i = 1:params.mesh.ngl# which basis

                    J = connijk[iel,i,j]

                    dψIJ_dx = dψ[i,k]*ψ[j,l]*dξdx[iel,k,l] + ψ[i,k]*dψ[j,l]*dηdx[iel,k,l]
                    dψIJ_dy = dψ[i,k]*ψ[j,l]*dξdy[iel,k,l] + ψ[i,k]*dψ[j,l]*dηdy[iel,k,l]

                    F_data[ip,2] = F_data[ip,2] .+ params.poisson[J,1].*dψIJ_dy
                    F_data[ip,1] = F_data[ip,1] .+ params.poisson[J,1].*dψIJ_dx
                
                end
            end

            count[ip] = count[ip] + 1

        end
    end

    F_data[:,2] = F_data[:,2]./count
    F_data[:,1] = F_data[:,1]./count

end

function poisson(params, connijk, x, y, ψ, RHS_laplacian_p, RHS_p)

    uind = unique(params.mesh.poin_in_bdy_edge)

    Laplacian(params, connijk, x, y, @view(params.uaux[:,:]))

    D_back = copy(@view(RHS_laplacian_p[uind,:,1]))
    b_back = copy(@view(RHS_p[uind,1]))

    boundary_poisson(params.neqs, uind, @view(RHS_laplacian_p[:,:,1]), @view(RHS_p[:,1]))    

    ψ[:,1] = compute_coef(@view(RHS_laplacian_p[:,:,:]), @view(RHS_p[:,:])) # a bit slow!

    RHS_laplacian_p[uind,:,1] .= D_back
    RHS_p[uind,1] .= b_back

end

function Laplacian(params, connijk, x, y, omega)

   resetRHS_pToZero_viscous!(params, params.SD)

    for ieq = 1:params.neqs
        for iel = 1:params.mesh.nelem
            compute_matirx_stiff(@view(params.rhs_laplacian_el_p[:,:,:,:]), 
                            ieq, 
                            params.ω, params.basis.ψ, params.basis.dψ, 
                            params.metrics.dξdx, params.metrics.dξdy, params.metrics.dηdx, params.metrics.dηdy, 
                            params.mesh.ngl, params.metrics.Je, 
                            connijk, iel)
            compute_matirx_rhs(@view(params.rhs_el_p[:,:,:]), 
                            ieq, 
                            params.ω, params.basis.ψ,
                            params.mesh.ngl, params.metrics.Je, 
                            connijk, iel, omega)
        end
    end

    DSS(@view(params.RHS_laplacian_p[:,:,:]), @view(params.rhs_laplacian_el_p[:,:,:,:]), 
        @view(params.RHS_p[:,:]), @view(params.rhs_el_p[:,:,:]), x, y, 
        connijk, params.neqs, params.mesh.nelem, params.mesh.ngl, 
        params.SD, params.AD)

end

function Laplacian2(params, connijk, x, y, omega)

    resetRHS_pToZero_viscous!(params, params.SD)

    for ieq = 1:params.neqs
        for iel = 1:params.mesh.nelem
            compute_matirx_stiff(@view(params.rhs_laplacian_el_p[:,:,:,:]), 
                            ieq, 
                            params.ω, params.basis.ψ, params.basis.dψ, 
                            params.metrics.dξdx, params.metrics.dξdy, params.metrics.dηdx, params.metrics.dηdy, 
                            params.mesh.ngl, params.metrics.Je, 
                            connijk, iel)
        end
    end

    DSS2(@view(params.RHS_laplacian_p[:,:,:]), @view(params.rhs_laplacian_el_p[:,:,:,:]),  
        connijk, params.neqs, params.mesh.nelem, params.mesh.ngl, 
        params.SD, params.AD)

end

function compute_boundary(params, Coef, u, t, nx, ny, x, y, connijk, dxdξ, dydξ)

    poin_in_bdy_edge = params.mesh.poin_in_bdy_edge
    nedges_bdy = params.mesh.nedges_bdy
    bdy_edge_in_elem = params.mesh.bdy_edge_in_elem
    bdy_edge_type = params.mesh.bdy_edge_type

    ngl = params.mesh.ngl
    neqs = params.neqs
    npoin = params.mesh.npoin

    ymax = params.mesh.ymax

    ω = params.ω

    if (t == 0)
        Laplacian2(params, connijk, x, y, @view(params.uaux[:,:]))
    end
    #K = @view(params.RHS_laplacian_p[:,:,1])

    velocity = zeros(nedges_bdy,ngl,2)

    resetSeg_pToZero_viscous!(params, params.SD)

    for iedge = 1:nedges_bdy
        for i=1:ngl

            ip = poin_in_bdy_edge[iedge,i]
            S = velocity[iedge, i, :]
            user_bc_dirichlet!(x[ip], y[ip], t, bdy_edge_type[iedge], S)
            velocity[iedge, i, :] = S

            I = poin_in_bdy_edge[iedge,i]
            iel  = bdy_edge_in_elem[iedge]

            if (y[I] == ymax)
                for l = 1:ngl, k = 1:ngl
                    ip = connijk[iel,k,l]

                    if (ip == I)
                        dψdn = -1*velocity[iedge, i, 1] + nx[iedge,i]*velocity[iedge, i, 2] # ???
                        J = sqrt((dxdξ[iel,k,l])^2 + (dydξ[iel,k,l]) ^2)
                        ωJ = ω[k]*J
                        @view(params.seg_p[I]) .+= ωJ * dψdn
                        break
                    end

                end
            end  

        end
    end

    w = params.Minv.*(@view(params.seg_p[:]) - @view(params.RHS_laplacian_p[:,:,1])*Coef)

    params.uaux[vec(poin_in_bdy_edge),1] = w[vec(poin_in_bdy_edge)]
    params.RHS[vec(poin_in_bdy_edge), 1] .= 0
    uaux2u!(u, @view(params.uaux[:,:]), neqs, npoin)    
    
end

# Exact quadrature

# Build stiff matrix
function compute_matirx_stiff(D_el, ieq, ω, ψ, dψ, dξdx, dξdy, dηdx, dηdy, ngl, Je, connijk, iel)

    for j = 1:ngl
        for i = 1:ngl # variable

            J = i + (j-1)*ngl

            for n = 1:ngl
                for m = 1:ngl # trail function

                    I = m + (n-1)*ngl

                    for l = 1:ngl
                        for k = 1:ngl # quadrature points on each element

                            ωkl  = ω[k]*ω[l]
                            Jkle = Je[iel, k, l]
                            index = connijk[iel,k,l]

                            dψij_dx = dψ[i,k]*ψ[j,l]*dξdx[iel,k,l] + ψ[i,k]*dψ[j,l]*dηdx[iel,k,l]
                            dψij_dy = dψ[i,k]*ψ[j,l]*dξdy[iel,k,l] + ψ[i,k]*dψ[j,l]*dηdy[iel,k,l]

                            dψmn_dx = dψ[m,k]*ψ[n,l]*dξdx[iel,k,l] + ψ[m,k]*dψ[n,l]*dηdx[iel,k,l]
                            dψmn_dy = dψ[m,k]*ψ[n,l]*dξdy[iel,k,l] + ψ[m,k]*dψ[n,l]*dηdy[iel,k,l]

                            D_el[iel,I,J,ieq] = D_el[iel,I,J,ieq] .+ ωkl*Jkle*(dψij_dx*dψmn_dx + dψij_dy*dψmn_dy)

                        end
                    end
                end
            end
        end
    end

end

# Build right hand side
function compute_matirx_rhs(b_el, ieq, ω, ψ, ngl, Je, connijk, iel, omega)

    for n = 1:ngl
        for m = 1:ngl # trail function
            I = m + (n-1)*ngl 
            for l = 1:ngl
                for k = 1:ngl # quadrature points on each element

                    ωkl  = ω[k]*ω[l]
                    Jkle = Je[iel, k, l]
                    index = connijk[iel,k,l]
                    ψmn = ψ[m,k]*ψ[n,l]
                    b_el[iel,I,ieq] = b_el[iel,I,ieq] .+ ωkl*Jkle*(omega[index])*ψmn
                end
            end
        end
    end
end

function compute_coef(D_global, b_global)

    prob = LinearProblem(sparse(D_global[:,:,1]), -b_global[:,1])  
    #sol = solve(prob, KLUFactorization())

    return solve(prob, KLUFactorization())
   
end

# My DSS
function DSS(D_global, D_el, b_global, b_el, x, y, connijk, neqs, nelem, ngl, ::NSD_2D, ::ContGal)

    for ieq = 1:neqs

        for iel = 1:nelem
            for j = 1:ngl,i = 1:ngl

                I = connijk[iel,i,j]
                b_global[I,ieq] = b_global[I,ieq] + b_el[iel,i+(j-1)*ngl,ieq]

                for n = 1:ngl,m = 1:ngl
                    J = connijk[iel,m,n]
                    D_global[I,J,ieq] = D_global[I,J,ieq] + D_el[iel,i+(j-1)*ngl,m+(n-1)*ngl,ieq]
                end

            end
        end

    end

end

function DSS2(D_global, D_el, connijk, neqs, nelem, ngl, ::NSD_2D, ::ContGal)

    for ieq = 1:neqs

        for iel = 1:nelem
            for j = 1:ngl,i = 1:ngl
                I = connijk[iel,i,j]
                for n = 1:ngl,m = 1:ngl
                    J = connijk[iel,m,n]
                    D_global[I,J,ieq] = D_global[I,J,ieq] + D_el[iel,i+(j-1)*ngl,m+(n-1)*ngl,ieq]
                end

            end
        end

    end

end


# apply boundary condition of stream function (Poisson Equation)
function boundary_poisson(neqs, uind, D_global, b_global)

    for ieq = 1:neqs

        D_global[uind, :, ieq] .= 0
        D_global[CartesianIndex.(uind, uind), ieq] .= 1
        b_global[uind, ieq] .= 0

    end
    
end
