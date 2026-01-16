using SparseArrays

function resetRHS_pToZero_viscous_VP!(params, SD)
    fill!(@view(params.RHS_laplacian_p[:,:,:]),  zero(params.T))
    fill!(@view(params.rhs_laplacian_el_p[:,:,:,:]), zero(params.T))
    fill!(@view(params.rhs_el_p[:,:,:]), zero(params.T))
    fill!(@view(params.RHS_p[:,:]),     zero(params.T))
end

# E (self-generated electric field)
function compute_∇V!(params, connijk, connijk_extra, x, v, F_data, SD::NSD_1D)

    ngl = params.mesh.ngl
    nelem = params.mesh.nelem
    npoin = params.mesh.npoin

    dψ = params.basis.dψ
    dξdx = params.metrics.dξdx

    poisson_VP(params, connijk, connijk_extra, x, v, 
            @view(params.poisson[:,:]), 
            @view(params.RHS_laplacian_p[:,:,:]), @view(params.RHS_p[:,:]), SD)

    # compute derivative of stream function
    fill!(F_data,zero(params.T))
    count = zeros(TFloat, npoin)

    for iel = 1:nelem
        for k = 1:ngl
            ip = connijk[iel,k,1]

            for i = 1:ngl# basis
                J = connijk[iel,i,1]

                dψdx = dψ[i,k]*dξdx[iel,k,1]
                F_data[ip] = F_data[ip] .+ params.poisson[J,1].*dψdx
            end

            count[ip] = count[ip] + 1
        end
    end

    F_data = F_data./count

end

# H (external electric field)
function compute_H!(params, connijk, connijk_extra, x, v, G_data, time, SD::NSD_1D)

    ngl = params.mesh.ngl
    nelem = params.mesh.nelem
    npoin = params.mesh.npoin

    dψ = params.basis.dψ
    dξdx = params.metrics.dξdx

    poisson_VP_external(params, connijk, connijk_extra, x, v, 
            @view(params.poisson[:,:]), 
            @view(params.RHS_laplacian_p[:,:,:]), @view(params.RHS_p[:,:]), SD, time)

    fill!(G_data,zero(params.T))
    count = zeros(TFloat, npoin)

    for iel = 1:nelem
        for k = 1:ngl
            ip = connijk[iel,k,1]

            for i = 1:ngl# basis
                J = connijk[iel,i,1]

                dψdx = dψ[i,k]*dξdx[iel,k,1]
                G_data[ip] = G_data[ip] .+ params.poisson[J,1].*dψdx
            end

            count[ip] = count[ip] + 1
        end
    end

    G_data = G_data./count

end

function compute_∇V!(params, connijk, connijk_extra, x, v, F_data, SD::NSD_2D)

    ngl = params.mesh.ngl
    nelem = params.mesh.nelem
    npoin = params.mesh.npoin

    dψ = params.basis.dψ
    ψ = params.basis.ψ
    dξdx = params.metrics.dξdx
    dξdy = params.metrics.dξdy
    dηdx = params.metrics.dηdx
    dηdy = params.metrics.dηdy

    poisson_VP(params, connijk, connijk_extra, x, v, 
            @view(params.poisson[:,:]), 
            @view(params.RHS_laplacian_p[:,:,:]), @view(params.RHS_p[:,:]), SD)

    # compute derivative of stream function
    fill!(F_data,zero(params.T))
    count = zeros(TFloat, npoin)

    for iel = 1:nelem
        for l = 1:ngl,k = 1:ngl
            ip = connijk[iel,k,l]

            for j = 1:ngl
                for i = 1:ngl# which basis

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

function compute_∇V!(params, connijk, connijk_extra, x, v, F_data, SD::NSD_3D)

    ngl = params.mesh.ngl
    nelem = params.mesh.nelem
    npoin = params.mesh.npoin

    dψ = params.basis.dψ
    ψ = params.basis.ψ
    dξdx = params.metrics.dξdx
    dξdy = params.metrics.dξdy
    dξdz = params.metrics.dξdz
    dηdx = params.metrics.dηdx
    dηdy = params.metrics.dηdy
    dηdz = params.metrics.dηdz
    dζdx = params.metrics.dζdx
    dζdy = params.metrics.dζdy
    dζdz = params.metrics.dζdz

    poisson_VP(params, connijk, connijk_extra, x, v, 
            @view(params.poisson[:,:]), 
            @view(params.RHS_laplacian_p[:,:,:]), @view(params.RHS_p[:,:]), SD)

    # compute derivative of stream function
    fill!(F_data,zero(params.T))
    count = zeros(TFloat, npoin)

    for iel = 1:nelem
        for p = 1:ngl, l = 1:ngl,k = 1:ngl
            ip = connijk[iel,k,l,p]

            for q = 1:ngl, j = 1:ngl, i = 1:ngl# which basis

                    J = connijk[iel,i,j,q]

                    dψIJ_dx = dψ[i,k]*ψ[j,l]*ψ[q,p]*dξdx[iel,k,l,p] + ψ[i,k]*dψ[j,l]*ψ[q,p]*dηdx[iel,k,l,p] + ψ[i,k]*ψ[j,l]*dψ[q,p]*dζdx[iel,k,l,p]
                    dψIJ_dy = dψ[i,k]*ψ[j,l]*ψ[q,p]*dξdy[iel,k,l,p] + ψ[i,k]*dψ[j,l]*ψ[q,p]*dηdy[iel,k,l,p] + ψ[i,k]*ψ[j,l]*dψ[q,p]*dζdy[iel,k,l,p]
                    dψIJ_dz = dψ[i,k]*ψ[j,l]*ψ[q,p]*dξdz[iel,k,l,p] + ψ[i,k]*dψ[j,l]*ψ[q,p]*dηdz[iel,k,l,p] + ψ[i,k]*ψ[j,l]*dψ[q,p]*dζdz[iel,k,l,p]

                    F_data[ip,3] = F_data[ip,3] .+ params.poisson[J,1].*dψIJ_dz
                    F_data[ip,2] = F_data[ip,2] .+ params.poisson[J,1].*dψIJ_dy
                    F_data[ip,1] = F_data[ip,1] .+ params.poisson[J,1].*dψIJ_dx
            end

            count[ip] = count[ip] + 1

        end
    end

    F_data[:,3] = F_data[:,3]./count
    F_data[:,2] = F_data[:,2]./count
    F_data[:,1] = F_data[:,1]./count

end

# E
function poisson_VP(params, connijk, connijk_extra, x, v, V, RHS_laplacian_p, RHS_p, SD) 

    omega = arrange_omega(params.mesh.nelem, params.mesh_extra.nelem, params.mesh.ngl, params.mesh_extra.ngl, params.mesh.npoin, params.mesh_extra.npoin, 
                         params.ω_extra, params.metrics_extra.Je, @view(params.uaux[:,1]),
                         connijk, connijk_extra, SD)

    # solve poisson equation: ΔV = ω

    Laplacian_VP(params, connijk, x, v, omega, SD)

    V[:,1] = compute_coef_VP(@view(RHS_laplacian_p[:,:,:]), @view(RHS_p[:,:]))
end

# H
function poisson_VP_external(params, connijk, connijk_extra, x, v, H, RHS_laplacian_p, RHS_p, SD, time) 

    nelem = params.mesh.nelem
    nelem_extra = params.mesh_extra.nelem

    npoin = params.mesh.npoin
    npoin_extra = params.mesh_extra.npoin

    ngl = params.mesh.ngl
    ngl_extra = params.mesh_extra.ngl

    epsilon = 0.005
    beta = 0.2
    v_avg = 2.4

    f = zeros(TFloat, npoin*npoin_extra)

    for iel = 1:nelem, i = 1:ngl
        ip = connijk[iel,i,1]

        for iel_extra = 1:nelem_extra, m = 1:ngl_extra
            ip_extra = connijk_extra[iel_extra,m,1]

            ind = (ip-1)*npoin_extra + ip_extra

            V = v[ip_extra]
            X = x[ip] - (time - 0.5*params.inputs[:Δt])*0.5*V

            mu = ( exp(-(V-v_avg)^2/2) + exp(-(V+v_avg)^2/2) )/(2*sqrt(2*pi))
            cs = cos(beta*X)

            f[ind] = epsilon*cs*mu
        end
    end

    omega = arrange_omega(params.mesh.nelem, params.mesh_extra.nelem, params.mesh.ngl, params.mesh_extra.ngl, params.mesh.npoin, params.mesh_extra.npoin, 
                         params.ω_extra, params.metrics_extra.Je, f,
                         connijk, connijk_extra, SD)

    # solve poisson equation: ΔV = ω

    Laplacian_VP(params, connijk, x, v, omega, SD)

    H[:,1] = compute_coef_VP(@view(RHS_laplacian_p[:,:,:]), @view(RHS_p[:,:]))
end

function arrange_omega(nelem, nelem_extra, ngl, ngl_extra, npoin, npoin_extra, ω_extra, Je_extra, f,
                       connijk, connijk_extra, SD::NSD_1D)

    rho = zeros(TFloat, npoin,1)
    count = zeros(TFloat, npoin)
    
   for iel = 1:nelem
        for i = 1:ngl
            ip = connijk[iel,i,1]

            for iel_extra = 1:nelem_extra
                for m = 1:ngl_extra
                    ip_extra = connijk_extra[iel_extra,m,1]

                    ind = (ip-1)*(npoin_extra) + ip_extra
                    ωJac_extra = ω_extra[m]*Je_extra[iel_extra,m,1]
                    rho[ip] = rho[ip] + ωJac_extra*f[ind]
                end
            end

            count[ip] = count[ip] + 1
        end
    end

    return  -rho./count # ΔV = - ρ
    # ρ = ∫ f dv
end

function arrange_omega(nelem, nelem_extra, ngl, ngl_extra, npoin, npoin_extra, ω_extra, Je_extra, f,
                       connijk, connijk_extra, SD::NSD_2D)

    rho = zeros(TFloat, npoin,1)
    count = zeros(TFloat, npoin)
    
   for iel = 1:nelem
        for i = 1:ngl
            for j = 1:ngl
                ip = connijk[iel,i,j]

                for iel_extra = 1:nelem_extra
                    for m = 1:ngl_extra
                        for n = 1:ngl_extra
                            ip_extra = connijk_extra[iel_extra,m,n]

                            ind = (ip-1)*(npoin_extra) + ip_extra
                            ωJac_extra = ω_extra[m]*ω_extra[n]*Je_extra[iel_extra,m,n]
                            rho[ip] = rho[ip] + ωJac_extra*f[ind]
                        end
                    end
                end
            
                count[ip] = count[ip] + 1
            end
        end
    end

    return  -rho./count # ΔV = - ρ
    # ρ = ∫ f dv
end

function arrange_omega(nelem, nelem_extra, ngl, ngl_extra, npoin, npoin_extra, ω_extra, Je_extra, f,
                       connijk, connijk_extra, SD::NSD_3D)

    rho = zeros(TFloat, npoin,1)
    count = zeros(TFloat, npoin)
    
    for iel = 1:nelem
        for k = 1:ngl
            for i = 1:ngl
                for j = 1:ngl
                    ip = connijk[iel,i,j,k]

                    for iel_extra = 1:nelem_extra
                        for l = 1:ngl_extra
                            for m = 1:ngl_extra
                                for n = 1:ngl_extra
                                    ip_extra = connijk_extra[iel_extra,m,n,l]

                                    ind = (ip-1)*(npoin_extra) + ip_extra
                                    ωJac_extra = ω_extra[m]*ω_extra[n]*ω_extra[l]*Je_extra[iel_extra,m,n,l]
                                    rho[ip] = rho[ip] + ωJac_extra*f[ind]
                                end
                            end
                        end
                    end
                
                    count[ip] = count[ip] + 1
                end
            end
        end
    end

    return  -rho./count # ΔV = - ρ
    # ρ = ∫ f dv
end


function Laplacian_VP(params, connijk, x, v, omega, SD::NSD_1D)

   resetRHS_pToZero_viscous_VP!(params, SD)

    for ieq = 1:params.neqs
        for iel = 1:params.mesh.nelem
            compute_matirx_stiff_VP(@view(params.rhs_laplacian_el_p[:,:,:,:]), 
                            ieq, 
                            params.ω, params.basis.ψ, params.basis.dψ, 
                            params.metrics.dξdx,
                            params.mesh.ngl, params.metrics.Je, 
                            connijk, iel, SD)
            compute_matirx_rhs_VP(@view(params.rhs_el_p[:,:,:]), 
                            ieq, 
                            params.ω, params.basis.ψ,
                            params.mesh.ngl, params.metrics.Je, 
                            connijk, iel, omega, SD)
        end
    end

    DSS_VP(@view(params.RHS_laplacian_p[:,:,:]), @view(params.rhs_laplacian_el_p[:,:,:,:]), 
        @view(params.RHS_p[:,:]), @view(params.rhs_el_p[:,:,:]), x, 
        connijk, params.neqs, params.mesh.nelem, params.mesh.ngl, SD, params.AD)

end

function Laplacian_VP(params, connijk, x, v, omega, SD::NSD_2D)

   resetRHS_pToZero_viscous_VP!(params, SD)

    for ieq = 1:params.neqs
        for iel = 1:params.mesh.nelem
            compute_matirx_stiff_VP(@view(params.rhs_laplacian_el_p[:,:,:,:]), 
                            ieq, 
                            params.ω, params.basis.ψ, params.basis.dψ, 
                            params.metrics.dξdx,params.metrics.dξdy,
                            params.metrics.dηdx,params.metrics.dηdy,
                            params.mesh.ngl, params.metrics.Je, 
                            connijk, iel, SD)
            compute_matirx_rhs_VP(@view(params.rhs_el_p[:,:,:]), 
                            ieq, 
                            params.ω, params.basis.ψ,
                            params.mesh.ngl, params.metrics.Je, 
                            connijk, iel, omega, SD)
        end
    end

    DSS_VP(@view(params.RHS_laplacian_p[:,:,:]), @view(params.rhs_laplacian_el_p[:,:,:,:]), 
        @view(params.RHS_p[:,:]), @view(params.rhs_el_p[:,:,:]), x, 
        connijk, params.neqs, params.mesh.nelem, params.mesh.ngl, SD, params.AD)

end

function Laplacian_VP(params, connijk, x, v, omega, SD::NSD_3D)

   resetRHS_pToZero_viscous_VP!(params, SD)

    for ieq = 1:params.neqs
        for iel = 1:params.mesh.nelem
            compute_matirx_stiff_VP(@view(params.rhs_laplacian_el_p[:,:,:,:]), 
                            ieq, 
                            params.ω, params.basis.ψ, params.basis.dψ, 
                            params.metrics.dξdx,params.metrics.dξdy,params.metrics.dξdz,
                            params.metrics.dηdx,params.metrics.dηdy,params.metrics.dηdz,
                            params.metrics.dζdx,params.metrics.dζdy,params.metrics.dζdz,
                            params.mesh.ngl, params.metrics.Je, 
                            connijk, iel, SD)
            compute_matirx_rhs_VP(@view(params.rhs_el_p[:,:,:]), 
                            ieq, 
                            params.ω, params.basis.ψ,
                            params.mesh.ngl, params.metrics.Je, 
                            connijk, iel, omega, SD)
        end
    end

    DSS_VP(@view(params.RHS_laplacian_p[:,:,:]), @view(params.rhs_laplacian_el_p[:,:,:,:]), 
        @view(params.RHS_p[:,:]), @view(params.rhs_el_p[:,:,:]), x, 
        connijk, params.neqs, params.mesh.nelem, params.mesh.ngl, SD, params.AD)

end

# Exact quadrature

# Build stiff matrix
function compute_matirx_stiff_VP(D_el, ieq, ω, ψ, dψ, dξdx, ngl, Je, connijk, iel, SD::NSD_1D)

        for i = 1:ngl # variable
            for m = 1:ngl # trail function
                for k = 1:ngl # quadrature points on each element

                    ωkl  = ω[k]
                    Jkle = Je[iel, k, 1]

                    dψij_dx = dψ[i,k]*dξdx[iel,k,1]
                    dψmn_dx = dψ[m,k]*dξdx[iel,k,1]

                    D_el[iel,i,m,ieq] = D_el[iel,i,m,ieq] .+ ωkl*Jkle*(dψij_dx*dψmn_dx)

                end                        
            end     
        end

end

function compute_matirx_stiff_VP(D_el, ieq, ω, ψ, dψ, dξdx, dξdy, dηdx, dηdy, ngl, Je, connijk, iel, SD::NSD_2D)

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

function compute_matirx_stiff_VP(D_el, ieq, ω, ψ, dψ, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, ngl, Je, connijk, iel, SD::NSD_3D)

    for p = 1:ngl,j = 1:ngl,i = 1:ngl # variable

            J = i + (j-1)*ngl + (p-1)*ngl^2

            for o = 1:ngl,n = 1:ngl,m = 1:ngl # trail function

                    I = m + (n-1)*ngl + (o-1)*ngl^2

                    for q = 1:ngl,l = 1:ngl,k = 1:ngl # quadrature points on each element

                            ωkl  = ω[k]*ω[l]*ω[q]
                            Jkle = Je[iel, k, l, q]
                            index = connijk[iel,k,l,q]

                            dψij_dx = dψ[i,k]*ψ[j,l]*ψ[p,q]*dξdx[iel,k,l,q] + ψ[i,k]*dψ[j,l]*ψ[p,q]*dηdx[iel,k,l,q] + ψ[i,k]*ψ[j,l]*dψ[p,q]*dζdx[iel,k,l,q]
                            dψij_dy = dψ[i,k]*ψ[j,l]*ψ[p,q]*dξdy[iel,k,l,q] + ψ[i,k]*dψ[j,l]*ψ[p,q]*dηdy[iel,k,l,q] + ψ[i,k]*ψ[j,l]*dψ[p,q]*dζdy[iel,k,l,q]
                            dψij_dz = dψ[i,k]*ψ[j,l]*ψ[p,q]*dξdz[iel,k,l,q] + ψ[i,k]*dψ[j,l]*ψ[p,q]*dηdz[iel,k,l,q] + ψ[i,k]*ψ[j,l]*dψ[p,q]*dζdz[iel,k,l,q]

                            dψmn_dx = dψ[m,k]*ψ[n,l]*ψ[o,q]*dξdx[iel,k,l,q] + ψ[m,k]*dψ[n,l]*ψ[o,q]*dηdx[iel,k,l,q] + ψ[m,k]*ψ[n,l]*dψ[o,q]*dζdx[iel,k,l,q]
                            dψmn_dy = dψ[m,k]*ψ[n,l]*ψ[o,q]*dξdy[iel,k,l,q] + ψ[m,k]*dψ[n,l]*ψ[o,q]*dηdy[iel,k,l,q] + ψ[m,k]*ψ[n,l]*dψ[o,q]*dζdy[iel,k,l,q]
                            dψmn_dz = dψ[m,k]*ψ[n,l]*ψ[o,q]*dξdz[iel,k,l,q] + ψ[m,k]*dψ[n,l]*ψ[o,q]*dηdz[iel,k,l,q] + ψ[m,k]*ψ[n,l]*dψ[o,q]*dζdz[iel,k,l,q]

                            D_el[iel,I,J,ieq] = D_el[iel,I,J,ieq] .+ ωkl*Jkle*(dψij_dx*dψmn_dx + dψij_dy*dψmn_dy + + dψij_dz*dψmn_dz)

                    end
            end
    end

end

# Build right hand side
function compute_matirx_rhs_VP(b_el, ieq, ω, ψ, ngl, Je, connijk, iel, omega, SD::NSD_1D)

        for m = 1:ngl # trail function
                for k = 1:ngl # quadrature points on each element

                    ωk  = ω[k]
                    Jke = Je[iel, k, 1]
                    index = connijk[iel,k,1]
                    ψm = ψ[m,k]
                    b_el[iel,m,ieq] = b_el[iel,m,ieq] .+ ωk*Jke*(omega[index])*ψm
                end
        end
end

function compute_matirx_rhs_VP(b_el, ieq, ω, ψ, ngl, Je, connijk, iel, omega, SD::NSD_2D)

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

function compute_matirx_rhs_VP(b_el, ieq, ω, ψ, ngl, Je, connijk, iel, omega, SD::NSD_3D)

    for o = 1:ngl,n = 1:ngl,m = 1:ngl # trail function
        I = m + (n-1)*ngl + (o-1)*ngl^2
        for p = 1:ngl,l = 1:ngl,k = 1:ngl # quadrature points on each element

            ωkl  = ω[k]*ω[l]*ω[p]
            Jkle = Je[iel, k, l, p]
            index = connijk[iel,k,l, p]
            ψmn = ψ[m,k]*ψ[n,l]*ψ[o,p]
            b_el[iel,I,ieq] = b_el[iel,I,ieq] .+ ωkl*Jkle*(omega[index])*ψmn
        end
    end
end


function compute_coef_VP(D_global, b_global)
    # Remark: D_global is just the stiff matrix ∫ψ_i'ψ_j', if you want to use it to solve Poisson equation, you have to add a minus to it.

    prob = LinearProblem(sparse(-D_global[:,:,1]), b_global[:,1])  
    #sol = solve(prob, KLUFactorization())

    return solve(prob, KLUFactorization())
   
end

# My DSS: 1D
function DSS_VP(D_global, D_el, b_global, b_el, x, connijk, neqs, nelem, ngl, ::NSD_1D, ::ContGal)

    for ieq = 1:neqs
        for iel = 1:nelem
            for i = 1:ngl

                I = connijk[iel,i,1]
                b_global[I,ieq] = b_global[I,ieq] + b_el[iel,i,ieq]

                for m = 1:ngl
                    J = connijk[iel,m,1]
                    D_global[I,J,ieq] = D_global[I,J,ieq] + D_el[iel,i,m,ieq]
                end
            end
        end
    end

    #boundary_poisson(neqs, ind_xmin, ind_xmax, D_global, b_global)

end

# My DSS: 2D
function DSS_VP(D_global, D_el, b_global, b_el, x, connijk, neqs, nelem, ngl, ::NSD_2D, ::ContGal)

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

# My DSS: 3D
function DSS_VP(D_global, D_el, b_global, b_el, x, connijk, neqs, nelem, ngl, ::NSD_3D, ::ContGal)

    for ieq = 1:neqs

        for iel = 1:nelem
            for k = 1:ngl,j = 1:ngl,i = 1:ngl

                I = connijk[iel,i,j,k]
                b_global[I,ieq] = b_global[I,ieq] + b_el[iel,i+(j-1)*ngl+(k-1)*ngl^2,ieq]

                for o = 1:ngl,n = 1:ngl,m = 1:ngl
                    J = connijk[iel,m,n,o]
                    D_global[I,J,ieq] = D_global[I,J,ieq] + D_el[iel,i+(j-1)*ngl+(k-1)*ngl^2,m+(n-1)*ngl+(o-1)*ngl^2,ieq]
                end

            end
        end

    end

end


# apply boundary condition of stream function (Poisson Equation)
function boundary_poisson(neqs, ind_xmin, ind_xmax, D_global, b_global)

    for ieq = 1:neqs

        D_global[ind_xmin, :, ieq] .= D_global[ind_xmin, :, ieq] + D_global[ind_xmax, :, ieq]
        mask = trues(size(D_global,1))
        mask[ind_xmax] .= false
        D_global = view(D_global, mask, :)

        D_global[:, ind_xmin, ieq] .= D_global[:, ind_xmin, ieq] + D_global[:, ind_xmax, ieq]
        mask = trues(size(D_global,2))
        mask[ind_xmax] .= false
        D_global = view(D_global, :, mask)

        b_global[ind_xmin, ieq] .= b_global[ind_xmin, ieq] + b_global[ind_xmax, ieq]
        mask = trues(size(b_global,1))
        mask[ind_xmax] .= false
        b_global = view(b_global, mask, :)

    end
    
end
