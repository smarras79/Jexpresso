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

function compute_∇ψ!(params, connijk, x, y, F_data)

    xmin = params.xmin; xmax = params.xmax; ymax = params.ymax

    ngl = params.mesh.ngl
    nelem = params.mesh.nelem
    neqs = params.neqs
    ω = params.ω
    ω1 = ω
    ψ = params.basis.ψ
    ψ1 = ψ
    dψ = params.basis.dψ
    dψ1 = dψ
    dξdx = params.metrics.dξdx
    dξdy = params.metrics.dξdy
    dηdx = params.metrics.dηdx
    dηdy = params.metrics.dηdy

     Coef = poisson(params, connijk, x, y)

    # compute derivative of stream function
    # 1. average method
    fill!(F_data,zeros(Float64))
    count = zeros(TFloat,maximum(connijk))

    for iel = 1:nelem
        for l = 1:ngl,k = 1:ngl
            ip = connijk[iel,k,l]

            for j = 1:params.mesh.ngl
                for i = 1:params.mesh.ngl# which basis

                    J = connijk[iel,i,j]

                    dψIJ_dx = dψ[i,k]*ψ1[j,l]*dξdx[iel,k,l] + ψ[i,k]*dψ1[j,l]*dηdx[iel,k,l]
                    dψIJ_dy = dψ[i,k]*ψ1[j,l]*dξdy[iel,k,l] + ψ[i,k]*dψ1[j,l]*dηdy[iel,k,l]

                    F_data[ip,2] = F_data[ip,2] .+ Coef[J].*dψIJ_dy
                    F_data[ip,1] = F_data[ip,1] .+ Coef[J].*dψIJ_dx
                
                end
            end

            count[ip] = count[ip] + 1

        end
    end

    F_data[:,2] = F_data[:,2]./count
    F_data[:,1] = F_data[:,1]./count


   #=  # 2. projection method
    derivative_y2 = zeros(TFloat, maximum(connijk))
    derivative_x2 = zeros(TFloat, maximum(connijk))
    Mel = zeros(TFloat, nelem, ngl*ngl, ngl*ngl)
    compute_mass(Mel, params)

    Mass = zeros(TFloat, maximum(connijk), maximum(connijk))
    DSS_mass(Mass, Mel, params)

    B_elx = zeros(TFloat, nelem, ngl*ngl)
    B_ely = zeros(TFloat, nelem, ngl*ngl)
    compute_rhs(B_elx, B_ely, Coef, params)

    Bx = zeros(TFloat, maximum(connijk))
    By = zeros(TFloat, maximum(connijk))
    DSS_rhs(Bx, B_elx, By, B_ely, params)

    derivative_y2 = Mass \ By
    derivative_x2 = Mass \ Bx =#  


end

function poisson(params, connijk, x, y)

    neqs = params.neqs

    D_global, b_global = Laplacian(params, connijk, x, y, params.uaux)

    Coef = zeros(TFloat, maximum(connijk), neqs)

    compute_coef(Coef, params, D_global, b_global)

    params.qp.qe = Coef

    return Coef

end

function Laplacian(params, connijk, x, y, omega)

    ngl = params.mesh.ngl
    nelem = params.mesh.nelem
    neqs = params.neqs
    ω = params.ω
    ω1 = ω
    ψ = params.basis.ψ
    ψ1 = ψ
    dψ = params.basis.dψ
    dψ1 = dψ
    dξdx = params.metrics.dξdx
    dξdy = params.metrics.dξdy
    dηdx = params.metrics.dηdx
    dηdy = params.metrics.dηdy

    D_el = zeros(TFloat, nelem, ngl*ngl, ngl*ngl, neqs)
    b_el = zeros(TFloat, nelem, ngl*ngl, neqs)

    for ieq = 1:params.neqs
        for iel = 1:params.mesh.nelem
            compute_matirx2(D_el,b_el,ieq,params,connijk,iel,x,y,omega)
        end
    end

    D_global = zeros(TFloat, maximum(connijk), maximum(connijk), params.neqs)
    b_global = zeros(TFloat, maximum(connijk), params.neqs)

    DSS(D_global, D_el, b_global, b_el, x, y, params, params.SD, params.AD)

    return D_global, b_global

end

function Laplacian2(params, connijk, x, y, omega)

    ngl = params.mesh.ngl
    nelem = params.mesh.nelem
    neqs = params.neqs
    ω = params.ω
    ω1 = ω
    ψ = params.basis.ψ
    ψ1 = ψ
    dψ = params.basis.dψ
    dψ1 = dψ
    dξdx = params.metrics.dξdx
    dξdy = params.metrics.dξdy
    dηdx = params.metrics.dηdx
    dηdy = params.metrics.dηdy

    D_el = zeros(TFloat, nelem, ngl*ngl, ngl*ngl, neqs)
    b_el = zeros(TFloat, nelem, ngl*ngl, neqs)

    for ieq = 1:params.neqs
        for iel = 1:params.mesh.nelem
            compute_matirx2(D_el,b_el,ieq,params,connijk,iel,x,y,omega)
        end
    end

    D_global = zeros(TFloat, maximum(connijk), maximum(connijk), params.neqs)
    b_global = zeros(TFloat, maximum(connijk), params.neqs)

    DSS2(D_global, D_el, b_global, b_el, x, y, params, params.SD, params.AD)

    return D_global, b_global

end

function compute_boundary(params, Coef, u)

    poin_in_bdy_edge = params.mesh.poin_in_bdy_edge
    nedges_bdy = params.mesh.nedges_bdy
    bdy_edge_in_elem = params.mesh.bdy_edge_in_elem

    nelem = params.mesh.nelem
    ngl = params.mesh.ngl
    neqs = params.neqs
    npoin = params.mesh.npoin

    connijk = params.mesh.connijk
    x = params.mesh.x
    y = params.mesh.y
    xmin = params.xmin
    xmax = params.xmax
    ymin = params.ymin
    ymax = params.ymax
    nx = params.metrics.nx
    ny = params.metrics.ny

    ω = params.ω
    ψ = params.basis.ψ
    ψ1 = ψ
    dψ = params.basis.dψ
    dψ1 = dψ

    dxdξ = params.metrics.dxdξ
    dydξ = params.metrics.dydξ
    dxdη = params.metrics.dxdη
    dydη = params.metrics.dydη

    IP = Vector{Float64}()
    El = Vector{Float64}()

    for iedge = 1:nedges_bdy 
        iel  = bdy_edge_in_elem[iedge]
        for k=1:ngl
            ip = poin_in_bdy_edge[iedge,k]
            push!(IP, ip)
            push!(El, iel)
        end
    end

    array = cat(IP', El', dims = 1)
    indices = unique(i -> array[1, i], 1:size(array, 2))
    array = array[:,indices]
    Nb = size(array)[2]
    index = Int.(array[1,:])

    D_global, b_global = Laplacian2(params, connijk, x, y, params.uaux)
    K = D_global[:,:,1]
    invM = params.Minv

    velocity = zeros(nedges_bdy,ngl,2)

    for iedge = 1:nedges_bdy
        for k=1:ngl
            user_bc_dirichlet!(x[ip], y[ip], t, bdy_edge_type[iedge], velocity[iedge,k,:])
        end
    end

    seg = zeros(TFloat, maximum(connijk))

    for iedge = 1:nedges_bdy
        for i = 1:ngl# basis

            I = poin_in_bdy[iedge,i]

            if (y[I] == ymax)
                for l = 1:ngl, k = 1:ngl
                    ip = connijk[iel,k,l]
                    
                    dψdn = ny[iedge,i]*velocity[iedge,i,1] + nx[iedge,i]*velocity[iedge,i,2]
                    J = sqrt((dxdξ[iel,k,l])^2 + (dydξ[iel,k,l]) ^2)
                    ωJ = ω[k]*J
                    seg[I] = seg[I] + ωJ * dψdn * ψ[i,k]*ψ1[j,l]
                end
            end           
        end
    end

    w = invM.*(seg - K*Coef)

    params.uaux[index,1] = w[index]
    params.RHS[index, 1] = zeros(TFloat, Nb)
    uaux2u!(u, params.uaux, neqs, npoin)    
end

function compute_oundary(params, Coef, u, index)

    connijk = params.mesh.connijk
    x = params.mesh.x
    y = params.mesh.y
    xmin = params.xmin
    xmax = params.xmax
    ymin = params.ymin
    ymax = params.ymax
    nelem = params.mesh.nelem
    ngl = params.mesh.ngl
    neqs = params.neqs
    npoin = params.mesh.npoin

    nbdy_edges = params.mesh.nedges_bdy
    bdy_edge_in_elem = params.mesh.bdy_edge_in_elem
    poin_in_bdy_edge = params.mesh.poin_in_bdy_edge

    ω = params.ω
    ω1 = ω
    ψ = params.basis.ψ
    ψ1 = ψ
    dψ = params.basis.dψ
    dψ1 = dψ

    dξdx = params.metrics.dξdx
    dξdy = params.metrics.dξdy
    dηdx = params.metrics.dηdx
    dηdy = params.metrics.dηdy

    dxdξ = params.metrics.dxdξ
    dydξ = params.metrics.dydξ
    dxdη = params.metrics.dxdη
    dydη = params.metrics.dydη

    D_global, b_global = Laplacian2(params, connijk, x, y, params.uaux)
    K = D_global[:,:,1]
    M = params.M
    invM = params.Minv

    velocity = 1

    seg = zeros(TFloat, maximum(connijk))

    for iel = 1:nelem
        for j = 1:ngl,i = 1:ngl # basis

            I = connijk[iel,i,j]

            if (y[I] == ymax)
                for l = 1:ngl, k = 1:ngl
                    ip = connijk[iel,k,l]

                    if (y[ip] == ymax)
                        dψdn = (-1*velocity)
                        J = sqrt((dxdξ[iel,k,l])^2 + (dydξ[iel,k,l]) ^2)
                        ωJ = ω[k]*J
                        seg[I] = seg[I] + ωJ*dψdn*ψ[i,k]*ψ1[j,l]
 
                    end

                end
            end           
        end
    end

    w = invM.*(seg - K*Coef)

    return w[index]

end

# Exact quadrature
function compute_matirx2(D_el,b_el,ieq,params,connijk,iel,x,y,omega)

    ω = params.ω
    ω1 = ω
    ψ = params.basis.ψ
    ψ1 = ψ
    dψ = params.basis.dψ
    dψ1 = dψ
    dξdx = params.metrics.dξdx
    dξdy = params.metrics.dξdy
    dηdx = params.metrics.dηdx
    dηdy = params.metrics.dηdy
#= 
    De = zeros(TFloat,params.mesh.ngl*params.mesh.ngl,params.mesh.ngl*params.mesh.ngl)
    be = zeros(TFloat,params.mesh.ngl*params.mesh.ngl,1)
 =#
    # Build stiff matrix
    for j = 1:params.mesh.ngl
        for i = 1:params.mesh.ngl # variable

            J = i + (j-1)*params.mesh.ngl

            for n = 1:params.mesh.ngl
                for m = 1:params.mesh.ngl # trail function

                    I = m + (n-1)*params.mesh.ngl

                    for l = 1:params.mesh.ngl
                        for k = 1:params.mesh.ngl # quadrature points on each element

                            ωkl  = ω[k]*ω1[l]
                            Jkle = params.metrics.Je[iel, k, l]
                            index = connijk[iel,k,l]

                            dψij_dx = dψ[i,k]*ψ1[j,l]*dξdx[iel,k,l] + ψ[i,k]*dψ1[j,l]*dηdx[iel,k,l]
                            dψij_dy = dψ[i,k]*ψ1[j,l]*dξdy[iel,k,l] + ψ[i,k]*dψ1[j,l]*dηdy[iel,k,l]

                            dψmn_dx = dψ[m,k]*ψ1[n,l]*dξdx[iel,k,l] + ψ[m,k]*dψ1[n,l]*dηdx[iel,k,l]
                            dψmn_dy = dψ[m,k]*ψ1[n,l]*dξdy[iel,k,l] + ψ[m,k]*dψ1[n,l]*dηdy[iel,k,l]

                            D_el[iel,I,J,ieq] = D_el[iel,I,J,ieq] .+ ωkl*Jkle*(dψij_dx*dψmn_dx + dψij_dy*dψmn_dy)

                        end
                    end
                end
            end
        end
    end

    # Build right hand side
            for n = 1:params.mesh.ngl
                for m = 1:params.mesh.ngl # trail function

                    I = m + (n-1)*params.mesh.ngl
                   
                    for l = 1:params.mesh.ngl
                        for k = 1:params.mesh.ngl # quadrature points on each element

                            ωkl  = ω[k]*ω1[l]
                            Jkle = params.metrics.Je[iel, k, l]
                            index = connijk[iel,k,l]
                            ψmn = ψ[m,k]*ψ1[n,l]

                            b_el[iel,I,ieq] = b_el[iel,I,ieq] .+ ωkl*Jkle*(omega[index])*ψmn

                        end
                    end
                end
            end

#= 
            # apply boundary condition for stream funtion
            for n = 1:params.mesh.ngl
                for m = 1:params.mesh.ngl

                    pos = connijk[iel,m,n]

                    I = m+(n-1)*params.mesh.ngl

                    if(x[pos] >= params.xmax || x[pos] <= params.xmin || y[pos] >= params.ymax || y[pos] <= params.ymin)                                                  
                        D_el[iel,I,:,ieq] = zeros(TFloat,1,params.mesh.ngl*params.mesh.ngl)
                        D_el[iel,I,I,ieq] = 1.0
                        be[I] = 0.0
                    end
                                                 
                end
            end
 =#
end


function compute_coef(Coef, params, D_global, b_global)

    neqs = params.neqs
    nelem = params.mesh.nelem
    ngl = params.mesh.ngl
    connijk = params.mesh.connijk

    for ieq = 1:neqs
        λ = 1e-9
        A = D_global[:,:,ieq] + λ*I
        b = -b_global[:,ieq]
        Coef[:,ieq] = A \ b
    end
   
end

# My DSS
function DSS(D_global, D_el, b_global, b_el, x, y, params, ::NSD_2D, ::ContGal)

    neqs = params.neqs
    nelem = params.mesh.nelem
    ngl = params.mesh.ngl
    connijk = params.mesh.connijk
    nx = params.metrics.nx
    ny = params.metrics.ny
    poin_in_bdy_edge = params.mesh.poin_in_bdy_edge
    nedges_bdy = params.mesh.nedges_bdy
    bdy_edge_in_elem = params.mesh.bdy_edge_in_elem

    xmin = params.xmin
    xmax = params.xmax
    ymin = params.ymin
    ymax = params.ymax
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
    dxdξ    = params.metrics.dxdξ
    connijk = params.mesh.connijk
        

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

        # apply boundary condition of stream function to the global stiff matrix
        for iel = 1:nelem
            for j = 1:ngl,i = 1:ngl # basis

                pos = connijk[iel,i,j]

                if (x[pos] >= xmax || x[pos] <= xmin || y[pos] <= ymin || y[pos] >= ymax)

                    D_global[pos,:,ieq] = zeros(TFloat,1,maximum(connijk))
                    D_global[pos,pos,ieq] = 1.0
                    b_global[pos,ieq] = 0.0

                end

               #=  if (y[pos] >= ymax)
                    iedge = findfirst(x -> x == iel, bdy_edge_in_elem)
                    Ny = zeros(TFloat, maximum(connijk))

                    for k = 1:ngl
                        ny_l = ny[iedge,k]
                        ip = poin_in_bdy_edge[iedge,k]
                        Ny[ip] = ny_l
                    end

                    for l = 1:ngl, k = 1:ngl
                        ip = connijk[iel,k,l]
                        if (y[ip] == ymax)
                            omega = ω[k]*ω1[l]
                            dψdn = Ny[ip]*(-1)
                            J = abs(dxdξ[iel,k,l])
                            b_global[pos,ieq] += omega*dψdn*ψ[i,k]*ψ1[j,l]*J
                        end
                    end
                end =#


               #=  if (y[pos] >= ymax)
                    u = 1
                    ip = find_ud(params, x[pos], y[pos], iel)
                    dy = abs(y[pos] - y[ip])
                    D_global[pos,:,ieq] = zeros(TFloat,1,maximum(connijk))
                    D_global[pos,pos,ieq] = 1.0
                    D_global[pos,ip,ieq] = -1.0
                    b_global[pos,ieq] = -dy*u
                end =#
                
            end
        end

    end

end

# DSS without boundary_condition
function DSS2(D_global, D_el, b_global, b_el, x, y, params, ::NSD_2D, ::ContGal)

    neqs = params.neqs
    nelem = params.mesh.nelem
    ngl = params.mesh.ngl
    connijk = params.mesh.connijk
    nx = params.metrics.nx
    ny = params.metrics.ny
    poin_in_bdy_edge = params.mesh.poin_in_bdy_edge
    nedges_bdy = params.mesh.nedges_bdy
    bdy_edge_in_elem = params.mesh.bdy_edge_in_elem

    xmin = params.xmin
    xmax = params.xmax
    ymin = params.ymin
    ymax = params.ymax
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
    dxdξ    = params.metrics.dxdξ
    connijk = params.mesh.connijk
        

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

function compute_rhs(B_elx, B_ely, Coef, params)
    ngl = params.mesh.ngl
    nelem = params.mesh.nelem
    neqs = params.neqs
    ω = params.ω
    ω1 = ω
    ψ = params.basis.ψ
    ψ1 = ψ
    dψ = params.basis.dψ
    dψ1 = dψ
    dξdx = params.metrics.dξdx
    dξdy = params.metrics.dξdy
    dηdx = params.metrics.dηdx
    dηdy = params.metrics.dηdy
    Je = params.metrics.Je
    connijk = params.mesh.connijk

    for iel = 1:nelem
        for n = 1:ngl,m = 1:ngl # basis

            tempy = 0.0;
            tempx = 0.0;

            for l = 1:ngl, k = 1:ngl # points

                ωkl  = ω[k]*ω1[l]
                Jkle = Je[iel, k, l]
                ψmn = ψ[m,k]*ψ1[n,l]

                for j = 1:params.mesh.ngl
                    for i = 1:params.mesh.ngl# which basis

                        J = connijk[iel,i,j]

                        dψIJ_dx = dψ[i,k]*ψ1[j,l]*dξdx[iel,k,l] + ψ[i,k]*dψ1[j,l]*dηdx[iel,k,l]
                        dψIJ_dy = dψ[i,k]*ψ1[j,l]*dξdy[iel,k,l] + ψ[i,k]*dψ1[j,l]*dηdy[iel,k,l]

                        tempy += Coef[J].*dψIJ_dy
                        tempx += Coef[J].*dψIJ_dx
                
                    end
                end

                B_elx[iel,m+(n-1)*ngl] += ωkl*Jkle*ψmn*tempx
                B_ely[iel,m+(n-1)*ngl] += ωkl*Jkle*ψmn*tempy

            end

        end
    end

end

function DSS_rhs(Bx, B_elx, By, B_ely, params)

    ngl = params.mesh.ngl
    nelem = params.mesh.nelem
    connijk = params.mesh.connijk

    for iel=1:nelem            
        for j = 1:ngl,i = 1:ngl
            I = connijk[iel,i,j]
            Bx[I] = Bx[I] + B_elx[iel,i+(j-1)*ngl]
            By[I] = By[I] + B_ely[iel,i+(j-1)*ngl]
        end
    end

end
