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
    count2 = Int64(params.number[1,1])
    print(time)
    print(" | ")
    print(maximum(velocity))
    print(" | ")
    print(count2)
    print(" | ")
    println(maximum(Coef))

    if (result == "velocity")
        data = velocity
    else
        data = Coef
    end

    # N = 2
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
        #if (count2>10 && count2<1200)

            #str = "result_$(lpad(count2, 5, '0'))"
            write_vtk_ref(SD, params.mesh, data, string(time), inputs[:output_dir]; nvar=length(params.qp.qe[1,:]), outvarsref=outvarsref)

        #end
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

function compute_boundary(params, Coef, u, t, nx, ny, x, y, connijk, dxdξ, dydξ, dxdη, dydη)

    poin_in_bdy_edge = params.mesh.poin_in_bdy_edge
    nedges_bdy = params.mesh.nedges_bdy
    bdy_edge_in_elem = params.mesh.bdy_edge_in_elem
    bdy_edge_type = params.mesh.bdy_edge_type

    ngl = params.mesh.ngl
    neqs = params.neqs
    npoin = params.mesh.npoin

    ymax = params.mesh.ymax
    ymin = params.mesh.ymin
    xmax = params.mesh.xmax
    xmin = params.mesh.xmin

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
            n_x = nx[iedge,i]
            n_y = ny[iedge,i]

            if (y[I] == ymax)
                n_y = 1
            elseif (y[I] == ymin)
                n_y = -1
            elseif (x[I] == xmax)
                n_x = 1
            elseif (x[I] == xmin)
                n_x = -1
            end

            for l = 1:ngl, k = 1:ngl
                ip = connijk[iel,k,l]
            #if (y[I] == ymax)
                if (ip == I)
                    dψdn = -n_y * velocity[iedge, i, 1] + n_x * velocity[iedge, i, 2]

                    if (y[I] == ymax || y[I] == ymin)
                        J = sqrt((dxdξ[iel,k,l])^2 + (dydξ[iel,k,l]) ^2)
                        ωJ = ω[k]*J
                    else
                        J = sqrt((dxdη[iel,k,l])^2 + (dydη[iel,k,l]) ^2)
                        ωJ = ω[l]*J
                    end

                    # J = sqrt((dxdξ[iel,k,l])^2 + (dydξ[iel,k,l]) ^2)        
                    # #J = abs(dxdξ[iel,k,l]*dydη[iel,k,l] - dxdη[iel,k,l]*dydξ[iel,k,l])
                    # ωJ = ω[k]*J

                    @view(params.seg_p[I]) .+= ωJ * dψdn
                    break
                end
            #end
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



# function draw_results(params, Coef, time, result)

#     SD      = params.SD
#     neqs    = params.neqs
#     ngl     = params.mesh.ngl
#     nelem   = params.mesh.nelem
#     ω       = params.ω
#     ω1      = ω
#     ψ       = params.basis.ψ
#     ψ1      = ψ
#     dψ      = params.basis.dψ
#     dψ1     = dψ
#     dξdx    = params.metrics.dξdx
#     dξdy    = params.metrics.dξdy
#     dηdx    = params.metrics.dηdx
#     dηdy    = params.metrics.dηdy
#     connijk = params.mesh.connijk 
#     x = params.mesh.x
#     y = params.mesh.y
#     xmin = params.xmin
#     xmax = params.xmax
#     ymin = params.ymin
#     ymax = params.ymax

#     inputs = params.inputs

#     velocity = zeros(TFloat,maximum(connijk))
#     velocity2 = zeros(TFloat,maximum(connijk))
#     count = zeros(TFloat,maximum(connijk))
#     for iel = 1:nelem
#         for j = 1:ngl, i = 1:ngl # points

#             ip = connijk[iel,i,j]

#             for n = 1:ngl, m = 1:ngl # basis functions

#                 J = connijk[iel,m,n]
#                 dψij_dy = dψ[m,i]*ψ1[n,j]*dξdy[iel,i,j] + ψ[m,i]*dψ1[n,j]*dηdy[iel,i,j]
#                 dψij_dx = dψ[m,i]*ψ1[n,j]*dξdx[iel,i,j] + ψ[m,i]*dψ1[n,j]*dηdx[iel,i,j]
#                 velocity[ip] = velocity[ip] - Coef[J]*dψij_dy
#                 velocity2[ip] = velocity2[ip] + Coef[J]*dψij_dx

#             end
#             count[ip] = count[ip] + 1
#         end
#     end
#     velocity = velocity./count
#     velocity2 = velocity2./count

#     outvarsref = Array{Union{Nothing, String}}(nothing, neqs)
#     qvars = ("ω")
#     for i = 1:length(outvarsref)
#         outvarsref[i] = string(qvars[i], "_ref")
#     end
    
#     print(time)
#     print(" | ")
#     print(maximum(velocity))
#     print(" | ")
#     println(maximum(Coef))

#     if (result == "velocity")
#         data = velocity
#     else
#         data = Coef
#     end

#     N = 2
#     res = zeros(TFloat, 3, maximum(connijk))

#     if (round(time,digits=3) == inputs[:tend])

#         for iel = 1:nelem
#             for j = 1:ngl, i = 1:ngl # points
#                 ip = connijk[iel,i,j]
#                 if (round(x[ip],digits=4) == (xmin+xmax)/2)
#                     res[1,ip] = y[ip]
#                     res[2,ip] = velocity[ip]
#                     res[3,ip] = velocity2[ip]
#                 end
#             end
#         end

#         ind = sortperm(reshape(res[1,:],(1,maximum(connijk))), dims = 2)
#         println("================================================")
#         println(res[2,ind])
#         println(res[3,ind])
#         println("================================================")

#         write_vtk_ref(SD, params.mesh, data, string(time), inputs[:output_dir]; nvar=length(params.qp.qe[1,:]), outvarsref=outvarsref)
#     end

# end

# function compute_∇ψ!(params, connijk, x, y, F_data)

#     xmin = params.xmin; xmax = params.xmax; ymax = params.ymax

#     ngl = params.mesh.ngl
#     nelem = params.mesh.nelem
#     neqs = params.neqs
#     ω = params.ω
#     ω1 = ω
#     ψ = params.basis.ψ
#     ψ1 = ψ
#     dψ = params.basis.dψ
#     dψ1 = dψ
#     dξdx = params.metrics.dξdx
#     dξdy = params.metrics.dξdy
#     dηdx = params.metrics.dηdx
#     dηdy = params.metrics.dηdy

#      Coef = poisson(params, connijk, x, y)

#     # compute derivative of stream function
#     # 1. average method
#     fill!(F_data,zero(params.T))
#     count = zeros(TFloat,maximum(connijk))

#     for iel = 1:nelem
#         for l = 1:ngl,k = 1:ngl
#             ip = connijk[iel,k,l]

#             for j = 1:params.mesh.ngl
#                 for i = 1:params.mesh.ngl# which basis

#                     J = connijk[iel,i,j]

#                     dψIJ_dx = dψ[i,k]*ψ1[j,l]*dξdx[iel,k,l] + ψ[i,k]*dψ1[j,l]*dηdx[iel,k,l]
#                     dψIJ_dy = dψ[i,k]*ψ1[j,l]*dξdy[iel,k,l] + ψ[i,k]*dψ1[j,l]*dηdy[iel,k,l]

#                     F_data[ip,2] = F_data[ip,2] .+ Coef[J].*dψIJ_dy
#                     F_data[ip,1] = F_data[ip,1] .+ Coef[J].*dψIJ_dx
                
#                 end
#             end

#             count[ip] = count[ip] + 1

#         end
#     end

#     F_data[:,2] = F_data[:,2]./count
#     F_data[:,1] = F_data[:,1]./count

# end

# function poisson(params, connijk, x, y)

#     neqs = params.neqs

#     D_global, b_global = Laplacian(params, connijk, x, y, params.uaux)

#     Coef = zeros(TFloat, maximum(connijk), neqs)

#     compute_coef(Coef, params, D_global, b_global)

#     params.qp.qe = Coef

#     return Coef

# end

# function Laplacian(params, connijk, x, y, omega)

#     ngl = params.mesh.ngl
#     nelem = params.mesh.nelem
#     neqs = params.neqs
#     ω = params.ω
#     ω1 = ω
#     ψ = params.basis.ψ
#     ψ1 = ψ
#     dψ = params.basis.dψ
#     dψ1 = dψ
#     dξdx = params.metrics.dξdx
#     dξdy = params.metrics.dξdy
#     dηdx = params.metrics.dηdx
#     dηdy = params.metrics.dηdy

#     D_el = zeros(TFloat, nelem, ngl*ngl, ngl*ngl, neqs)
#     b_el = zeros(TFloat, nelem, ngl*ngl, neqs)

#     for ieq = 1:params.neqs
#         for iel = 1:params.mesh.nelem
#             compute_matirx2(D_el,b_el,ieq,params,connijk,iel,x,y,omega)
#         end
#     end

#     D_global = zeros(TFloat, maximum(connijk), maximum(connijk), params.neqs)
#     b_global = zeros(TFloat, maximum(connijk), params.neqs)

#     DSS(D_global, D_el, b_global, b_el, x, y, params, params.SD, params.AD)

#     return D_global, b_global

# end

# function Laplacian2(params, connijk, x, y, omega)

#     ngl = params.mesh.ngl
#     nelem = params.mesh.nelem
#     neqs = params.neqs
#     ω = params.ω
#     ω1 = ω
#     ψ = params.basis.ψ
#     ψ1 = ψ
#     dψ = params.basis.dψ
#     dψ1 = dψ
#     dξdx = params.metrics.dξdx
#     dξdy = params.metrics.dξdy
#     dηdx = params.metrics.dηdx
#     dηdy = params.metrics.dηdy

#     D_el = zeros(TFloat, nelem, ngl*ngl, ngl*ngl, neqs)
#     b_el = zeros(TFloat, nelem, ngl*ngl, neqs)

#     for ieq = 1:params.neqs
#         for iel = 1:params.mesh.nelem
#             compute_matirx2(D_el,b_el,ieq,params,connijk,iel,x,y,omega)
#         end
#     end

#     D_global = zeros(TFloat, maximum(connijk), maximum(connijk), params.neqs)
#     b_global = zeros(TFloat, maximum(connijk), params.neqs)

#     DSS2(D_global, D_el, b_global, b_el, x, y, params, params.SD, params.AD)

#     return D_global, b_global

# end

# function compute_boundary(params, Coef, u, t)

#     poin_in_bdy_edge = params.mesh.poin_in_bdy_edge
#     nedges_bdy = params.mesh.nedges_bdy
#     bdy_edge_in_elem = params.mesh.bdy_edge_in_elem
#     bdy_edge_type = params.mesh.bdy_edge_type

#     nelem = params.mesh.nelem
#     ngl = params.mesh.ngl
#     neqs = params.neqs
#     npoin = params.mesh.npoin

#     connijk = params.mesh.connijk
#     x = params.mesh.x
#     y = params.mesh.y
#     xmin = params.xmin
#     xmax = params.xmax
#     ymin = params.ymin
#     ymax = params.ymax
#     nx = params.metrics.nx
#     ny = params.metrics.ny

#     ω = params.ω
#     ψ = params.basis.ψ
#     ψ1 = ψ
#     dψ = params.basis.dψ
#     dψ1 = dψ

#     dxdξ = params.metrics.dxdξ
#     dydξ = params.metrics.dydξ
#     dxdη = params.metrics.dxdη
#     dydη = params.metrics.dydη

#     IP = Vector{Float64}()
#     El = Vector{Float64}()

#     for iedge = 1:nedges_bdy 
#         iel  = bdy_edge_in_elem[iedge]
#         for k=1:ngl
#             ip = poin_in_bdy_edge[iedge,k]
#             push!(IP, ip)
#             push!(El, iel)
#         end
#     end

#     array = cat(IP', El', dims = 1)
#     indices = unique(i -> array[1, i], 1:size(array, 2))
#     array = array[:,indices]
#     Nb = size(array)[2]
#     index = Int.(array[1,:])

#     D_global, b_global = Laplacian2(params, connijk, x, y, params.uaux)
#     K = D_global[:,:,1]
#     invM = params.Minv

#     velocity = zeros(nedges_bdy,ngl,2)

#     for iedge = 1:nedges_bdy
#         for k=1:ngl
#             ip = poin_in_bdy_edge[iedge,k]
#             S = velocity[iedge,k,:]
#             user_bc_dirichlet!(x[ip], y[ip], t, bdy_edge_type[iedge], S)
#             velocity[iedge,k,:] = S
#         end
#     end

#     seg = zeros(TFloat, maximum(connijk))

#     for iedge = 1:nedges_bdy
#         for i = 1:ngl# basis

#             I = poin_in_bdy_edge[iedge,i]
#             iel = bdy_edge_in_elem[iedge]

#             if (y[I] == ymax)
#                 for l = 1:ngl, k = 1:ngl
#                     ip = connijk[iel,k,l]
#                     if (ip == I)
                        
                    
#                         dψdn = -velocity[iedge,k,1]
#                         J = sqrt((dxdξ[iel,k,l])^2 + (dydξ[iel,k,l]) ^2)
#                         ωJ = ω[k]*J
#                         seg[I] = seg[I] + ωJ * dψdn
#                         break
#                     end
#                 end
#             end           
#         end
#     end

#     w = invM.*(seg - K*Coef)

#     params.uaux[index,1] = w[index]
#     params.RHS[index, 1] = zeros(TFloat, Nb)
#     uaux2u!(u, params.uaux, neqs, npoin)    
# end

# # Exact quadrature
# function compute_matirx2(D_el,b_el,ieq,params,connijk,iel,x,y,omega)

#     ω = params.ω
#     ω1 = ω
#     ψ = params.basis.ψ
#     ψ1 = ψ
#     dψ = params.basis.dψ
#     dψ1 = dψ
#     dξdx = params.metrics.dξdx
#     dξdy = params.metrics.dξdy
#     dηdx = params.metrics.dηdx
#     dηdy = params.metrics.dηdy
# #= 
#     De = zeros(TFloat,params.mesh.ngl*params.mesh.ngl,params.mesh.ngl*params.mesh.ngl)
#     be = zeros(TFloat,params.mesh.ngl*params.mesh.ngl,1)
#  =#
#     # Build stiff matrix
#     for j = 1:params.mesh.ngl
#         for i = 1:params.mesh.ngl # variable

#             J = i + (j-1)*params.mesh.ngl

#             for n = 1:params.mesh.ngl
#                 for m = 1:params.mesh.ngl # trail function

#                     I = m + (n-1)*params.mesh.ngl

#                     for l = 1:params.mesh.ngl
#                         for k = 1:params.mesh.ngl # quadrature points on each element

#                             ωkl  = ω[k]*ω1[l]
#                             Jkle = params.metrics.Je[iel, k, l]
#                             index = connijk[iel,k,l]

#                             dψij_dx = dψ[i,k]*ψ1[j,l]*dξdx[iel,k,l] + ψ[i,k]*dψ1[j,l]*dηdx[iel,k,l]
#                             dψij_dy = dψ[i,k]*ψ1[j,l]*dξdy[iel,k,l] + ψ[i,k]*dψ1[j,l]*dηdy[iel,k,l]

#                             dψmn_dx = dψ[m,k]*ψ1[n,l]*dξdx[iel,k,l] + ψ[m,k]*dψ1[n,l]*dηdx[iel,k,l]
#                             dψmn_dy = dψ[m,k]*ψ1[n,l]*dξdy[iel,k,l] + ψ[m,k]*dψ1[n,l]*dηdy[iel,k,l]

#                             D_el[iel,I,J,ieq] = D_el[iel,I,J,ieq] .+ ωkl*Jkle*(dψij_dx*dψmn_dx + dψij_dy*dψmn_dy)

#                         end
#                     end
#                 end
#             end
#         end
#     end

#     # Build right hand side
#             for n = 1:params.mesh.ngl
#                 for m = 1:params.mesh.ngl # trail function

#                     I = m + (n-1)*params.mesh.ngl
                   
#                     for l = 1:params.mesh.ngl
#                         for k = 1:params.mesh.ngl # quadrature points on each element

#                             ωkl  = ω[k]*ω1[l]
#                             Jkle = params.metrics.Je[iel, k, l]
#                             index = connijk[iel,k,l]
#                             ψmn = ψ[m,k]*ψ1[n,l]

#                             b_el[iel,I,ieq] = b_el[iel,I,ieq] .+ ωkl*Jkle*(omega[index])*ψmn

#                         end
#                     end
#                 end
#             end

# #= 
#             # apply boundary condition for stream funtion
#             for n = 1:params.mesh.ngl
#                 for m = 1:params.mesh.ngl

#                     pos = connijk[iel,m,n]

#                     I = m+(n-1)*params.mesh.ngl

#                     if(x[pos] >= params.xmax || x[pos] <= params.xmin || y[pos] >= params.ymax || y[pos] <= params.ymin)                                                  
#                         D_el[iel,I,:,ieq] = zeros(TFloat,1,params.mesh.ngl*params.mesh.ngl)
#                         D_el[iel,I,I,ieq] = 1.0
#                         be[I] = 0.0
#                     end
                                                 
#                 end
#             end
#  =#
# end


# function compute_coef(Coef, params, D_global, b_global)

#     neqs = params.neqs
#     nelem = params.mesh.nelem
#     ngl = params.mesh.ngl
#     connijk = params.mesh.connijk

#     for ieq = 1:neqs
#         λ = 1e-9
#         A = D_global[:,:,ieq] + λ*I
#         b = -b_global[:,ieq]
#         Coef[:,ieq] = A \ b
#     end
   
# end

# # My DSS
# function DSS(D_global, D_el, b_global, b_el, x, y, params, ::NSD_2D, ::ContGal)

#     neqs = params.neqs
#     nelem = params.mesh.nelem
#     ngl = params.mesh.ngl
#     connijk = params.mesh.connijk
#     nx = params.metrics.nx
#     ny = params.metrics.ny
#     poin_in_bdy_edge = params.mesh.poin_in_bdy_edge
#     nedges_bdy = params.mesh.nedges_bdy
#     bdy_edge_in_elem = params.mesh.bdy_edge_in_elem

#     xmin = params.xmin
#     xmax = params.xmax
#     ymin = params.ymin
#     ymax = params.ymax
#     ω       = params.ω
#     ω1      = ω
#     ψ       = params.basis.ψ
#     ψ1      = ψ
#     dψ      = params.basis.dψ
#     dψ1     = dψ
#     dξdx    = params.metrics.dξdx
#     dξdy    = params.metrics.dξdy
#     dηdx    = params.metrics.dηdx
#     dηdy    = params.metrics.dηdy
#     dxdξ    = params.metrics.dxdξ
#     connijk = params.mesh.connijk
        

#     for ieq = 1:neqs

#         for iel = 1:nelem
#             for j = 1:ngl,i = 1:ngl

#                 I = connijk[iel,i,j]
#                 b_global[I,ieq] = b_global[I,ieq] + b_el[iel,i+(j-1)*ngl,ieq]

#                 for n = 1:ngl,m = 1:ngl
#                     J = connijk[iel,m,n]
#                     D_global[I,J,ieq] = D_global[I,J,ieq] + D_el[iel,i+(j-1)*ngl,m+(n-1)*ngl,ieq]
#                 end

#             end
#         end

#         # apply boundary condition of stream function to the global stiff matrix
#         for iel = 1:nelem
#             for j = 1:ngl,i = 1:ngl # basis

#                 pos = connijk[iel,i,j]

#                 if (x[pos] >= xmax || x[pos] <= xmin || y[pos] <= ymin || y[pos] >= ymax)

#                     D_global[pos,:,ieq] = zeros(TFloat,1,maximum(connijk))
#                     D_global[pos,pos,ieq] = 1.0
#                     b_global[pos,ieq] = 0.0

#                 end

#                #=  if (y[pos] >= ymax)
#                     iedge = findfirst(x -> x == iel, bdy_edge_in_elem)
#                     Ny = zeros(TFloat, maximum(connijk))

#                     for k = 1:ngl
#                         ny_l = ny[iedge,k]
#                         ip = poin_in_bdy_edge[iedge,k]
#                         Ny[ip] = ny_l
#                     end

#                     for l = 1:ngl, k = 1:ngl
#                         ip = connijk[iel,k,l]
#                         if (y[ip] == ymax)
#                             omega = ω[k]*ω1[l]
#                             dψdn = Ny[ip]*(-1)
#                             J = abs(dxdξ[iel,k,l])
#                             b_global[pos,ieq] += omega*dψdn*ψ[i,k]*ψ1[j,l]*J
#                         end
#                     end
#                 end =#


#                #=  if (y[pos] >= ymax)
#                     u = 1
#                     ip = find_ud(params, x[pos], y[pos], iel)
#                     dy = abs(y[pos] - y[ip])
#                     D_global[pos,:,ieq] = zeros(TFloat,1,maximum(connijk))
#                     D_global[pos,pos,ieq] = 1.0
#                     D_global[pos,ip,ieq] = -1.0
#                     b_global[pos,ieq] = -dy*u
#                 end =#
                
#             end
#         end

#     end

# end

# # DSS without boundary_condition
# function DSS2(D_global, D_el, b_global, b_el, x, y, params, ::NSD_2D, ::ContGal)

#     neqs = params.neqs
#     nelem = params.mesh.nelem
#     ngl = params.mesh.ngl
#     connijk = params.mesh.connijk
#     nx = params.metrics.nx
#     ny = params.metrics.ny
#     poin_in_bdy_edge = params.mesh.poin_in_bdy_edge
#     nedges_bdy = params.mesh.nedges_bdy
#     bdy_edge_in_elem = params.mesh.bdy_edge_in_elem

#     xmin = params.xmin
#     xmax = params.xmax
#     ymin = params.ymin
#     ymax = params.ymax
#     ω       = params.ω
#     ω1      = ω
#     ψ       = params.basis.ψ
#     ψ1      = ψ
#     dψ      = params.basis.dψ
#     dψ1     = dψ
#     dξdx    = params.metrics.dξdx
#     dξdy    = params.metrics.dξdy
#     dηdx    = params.metrics.dηdx
#     dηdy    = params.metrics.dηdy
#     dxdξ    = params.metrics.dxdξ
#     connijk = params.mesh.connijk
        

#     for ieq = 1:neqs

#         for iel = 1:nelem
#             for j = 1:ngl,i = 1:ngl

#                 I = connijk[iel,i,j]
#                 b_global[I,ieq] = b_global[I,ieq] + b_el[iel,i+(j-1)*ngl,ieq]

#                 for n = 1:ngl,m = 1:ngl
#                     J = connijk[iel,m,n]
#                     D_global[I,J,ieq] = D_global[I,J,ieq] + D_el[iel,i+(j-1)*ngl,m+(n-1)*ngl,ieq]
#                 end

#             end
#         end
#     end

# end

# function compute_rhs(B_elx, B_ely, Coef, params)
#     ngl = params.mesh.ngl
#     nelem = params.mesh.nelem
#     neqs = params.neqs
#     ω = params.ω
#     ω1 = ω
#     ψ = params.basis.ψ
#     ψ1 = ψ
#     dψ = params.basis.dψ
#     dψ1 = dψ
#     dξdx = params.metrics.dξdx
#     dξdy = params.metrics.dξdy
#     dηdx = params.metrics.dηdx
#     dηdy = params.metrics.dηdy
#     Je = params.metrics.Je
#     connijk = params.mesh.connijk

#     for iel = 1:nelem
#         for n = 1:ngl,m = 1:ngl # basis

#             tempy = 0.0;
#             tempx = 0.0;

#             for l = 1:ngl, k = 1:ngl # points

#                 ωkl  = ω[k]*ω1[l]
#                 Jkle = Je[iel, k, l]
#                 ψmn = ψ[m,k]*ψ1[n,l]

#                 for j = 1:params.mesh.ngl
#                     for i = 1:params.mesh.ngl# which basis

#                         J = connijk[iel,i,j]

#                         dψIJ_dx = dψ[i,k]*ψ1[j,l]*dξdx[iel,k,l] + ψ[i,k]*dψ1[j,l]*dηdx[iel,k,l]
#                         dψIJ_dy = dψ[i,k]*ψ1[j,l]*dξdy[iel,k,l] + ψ[i,k]*dψ1[j,l]*dηdy[iel,k,l]

#                         tempy += Coef[J].*dψIJ_dy
#                         tempx += Coef[J].*dψIJ_dx
                
#                     end
#                 end

#                 B_elx[iel,m+(n-1)*ngl] += ωkl*Jkle*ψmn*tempx
#                 B_ely[iel,m+(n-1)*ngl] += ωkl*Jkle*ψmn*tempy

#             end

#         end
#     end

# end

# function DSS_rhs(Bx, B_elx, By, B_ely, params)

#     ngl = params.mesh.ngl
#     nelem = params.mesh.nelem
#     connijk = params.mesh.connijk

#     for iel=1:nelem            
#         for j = 1:ngl,i = 1:ngl
#             I = connijk[iel,i,j]
#             Bx[I] = Bx[I] + B_elx[iel,i+(j-1)*ngl]
#             By[I] = By[I] + B_ely[iel,i+(j-1)*ngl]
#         end
#     end

# end