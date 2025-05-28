function build_radiative_transfer_problem(mesh, inputs, neqs, ngl, dψ, ψ, ω, Je, dξdx, dηdy, extra_mesh::St_extra_mesh, QT::Inexact, SD::NSD_2D, AD::ContGal)

    npoin = mesh.npoin
    nelem = mesh.nelem
    #LHS_ang_spat = Array{Array,3}(undef, Int64(nelem), Int64(ngl^2), Int64(ngl^2))
    #Mass_ang_spat = Array{Array,3}(undef, Int64(nelem), Int64(ngl^2), Int64(ngl^2))
    if (inputs[:adaptive_extra_meshes])
        LHS_ang_spat = Array{Array,3}(undef, Int64(nelem), Int64(ngl^2), Int64(ngl^2))
        Mass_ang_spat = Array{Array,3}(undef, Int64(nelem), Int64(ngl^2), Int64(ngl^2))
        npoin_ang_total = 0
        for iel=1:nelem
            npoin_ang_total += extra_mesh[iel].extra_npoin*ngl*ngl
            nelem_ext = extra_mesh[iel].extra_nelem
            ngl_ext = extra_mesh[iel].extra_nop[1]+1
            κ = zeros(TFloat, ngl, ngl)            
            σ = zeros(TFloat, ngl, ngl)
            Φ = zeros(TFloat, extra_mesh[iel].extra_npoin, extra_mesh[iel].extra_npoin)
            E_rad_el = zeros(TFloat, nelem, ngl^2, ngl^2, nelem_ext, extra_mesh[iel].extra_nop[1]+1, extra_mesh[iel].extra_nop[1]+1)
            P_rad_el = zeros(TFloat, nelem, ngl^2, ngl^2, nelem_ext, extra_mesh[iel].extra_nop[1]+1, extra_mesh[iel].extra_nop[1]+1)
            S_rad_el = zeros(TFloat, nelem, ngl^2, ngl^2, nelem_ext, extra_mesh[iel].extra_nop[1]+1, extra_mesh[iel].extra_nop[1]+1)
            user_extinction!(κ, mesh.x, mesh.y, mesh.connijk, iel, extra_mesh, ngl)
            user_scattering!(σ, Φ, mesh.x, mesh.y, mesh.connijk, iel, extra_mesh, ngl)
            _rad_element_extinction_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, E_rad_el, iel, κ, extra_mesh, QT, SD, AD)
            _rad_element_internal_propagation_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, dξdx, dηdy, P_rad_el, iel, extra_mesh, QT, SD, AD)
            _rad_element_scattering_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, S_rad_el, iel, σ, Φ, extra_mesh, QT, SD, AD)
            for i=1:ngl
                for j =1:ngl
                    LHS_ang_spat[iel,i,j] = zeros(TFloat, extra_mesh[iel].extra_npoin, extra_mesh[iel].extra_npoin)
                end
            end
        end
        LHS_el_ang = E_rad_el + P_rad_el - S_rad_el 
        for iel=1:nelem
            DSS_angular!(LHS_el_ang, LHS_ang_spat, nelem, extra_mesh, ngl, mesh.connijk, iel, SD)
        end

    else
        LHS_ang_spat = zeros(Float64, Int64(nelem), Int64(ngl^2), Int64(ngl^2), extra_mesh.extra_npoin, extra_mesh.extra_npoin)
        Mass_ang_spat = zeros(Float64, Int64(nelem), Int64(ngl^2), Int64(ngl^2), extra_mesh.extra_npoin, extra_mesh.extra_npoin)
        npoin_ang_total = npoin*extra_mesh.extra_npoin
        nelem_ext = extra_mesh.extra_nelem
        E_rad_el = zeros(TFloat, nelem, ngl^2, ngl^2, nelem_ext, extra_mesh.extra_nop[1]+1, extra_mesh.extra_nop[1]+1)
        P_rad_el = zeros(TFloat, nelem, ngl^2, ngl^2, nelem_ext, extra_mesh.extra_nop[1]+1, extra_mesh.extra_nop[1]+1)
        S_rad_el = zeros(TFloat, nelem, ngl^2, ngl^2, nelem_ext, extra_mesh.extra_nop[1]+1, extra_mesh.extra_nop[1]+1)
        M_rad_el = zeros(TFloat, nelem, ngl^2, ngl^2, nelem_ext, extra_mesh.extra_nop[1]+1, extra_mesh.extra_nop[1]+1)
        κ = zeros(TFloat, ngl, ngl)
        σ = zeros(TFloat, ngl, ngl)
        Φ = zeros(TFloat, extra_mesh.extra_npoin, extra_mesh.extra_npoin)
        for iel=1:nelem
            fill!(κ, zero(Float64))
            fill!(σ, zero(Float64))
            fill!(Φ, zero(Float64))
            user_extinction!(κ, mesh.x, mesh.y, mesh.connijk, iel, extra_mesh, ngl)
            user_scattering!(σ, Φ, mesh.x, mesh.y, mesh.connijk, iel, extra_mesh, ngl)
            _rad_element_extinction_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, E_rad_el, iel, κ, extra_mesh.extra_nop, extra_mesh.extra_nelem, extra_mesh.extra_metrics.Je, extra_mesh.ωθ, 
                                                   extra_mesh.ψ, QT, SD, AD)
            _rad_element_internal_propagation_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, dξdx, dηdy, P_rad_el, iel, extra_mesh.extra_nop, extra_mesh.extra_nelem, extra_mesh.extra_metrics.Je, extra_mesh.ωθ,
                                                   extra_mesh.ψ, extra_mesh.extra_connijk, extra_mesh.extra_coords, QT, SD, AD)
            _rad_element_scattering_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, S_rad_el, iel, σ, Φ, extra_mesh.extra_nop, extra_mesh.extra_nelem, extra_mesh.extra_metrics.Je, extra_mesh.ωθ,
                                                   extra_mesh.ψ, extra_mesh.extra_connijk, QT, SD, AD)
            _rad_element_mass_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, M_rad_el, iel, extra_mesh.extra_nop, extra_mesh.extra_nelem, extra_mesh.extra_metrics.Je, extra_mesh.ωθ,
                                                   extra_mesh.ψ, QT, SD, AD)
        end
        LHS_el_ang = E_rad_el + P_rad_el - S_rad_el
        for iel=1:nelem
            DSS_angular!(LHS_el_ang, LHS_ang_spat, nelem, extra_mesh.extra_nop, extra_mesh.extra_nelem, extra_mesh.extra_connijk, ngl, mesh.connijk, iel, SD)
            DSS_angular!(M_rad_el, Mass_ang_spat, nelem, extra_mesh.extra_nop, extra_mesh.extra_nelem, extra_mesh.extra_connijk, ngl, mesh.connijk, iel, SD)
        end
        LHS_sp = zeros(TFloat, npoin, npoin, extra_mesh.extra_npoin, extra_mesh.extra_npoin)
        M_sp = zeros(TFloat, npoin, npoin, extra_mesh.extra_npoin, extra_mesh.extra_npoin)
        DSS_spatial!(LHS_sp, LHS_ang_spat, nelem, extra_mesh.extra_npoin, ngl, mesh.connijk, inputs, SD)
        DSS_spatial!(M_sp, Mass_ang_spat, nelem, extra_mesh.extra_npoin, ngl, mesh.connijk, inputs, SD)
    end
    LHS = zeros(TFloat, npoin_ang_total, npoin_ang_total)
    M = zeros(TFloat, npoin_ang_total, npoin_ang_total)
    Map_rad_global!(LHS_sp, LHS, npoin, extra_mesh.extra_npoin, nelem, mesh.connijk, ngl, inputs, SD)
    @info "assembled LHS"
    Map_rad_global!(M_sp, M, npoin, extra_mesh.extra_npoin, nelem, mesh.connijk, ngl, inputs, SD)
    @info "assembled Mass matrix"
    M_inv = zeros(TFloat, npoin_ang_total, npoin_ang_total)
    M_inv = inv(M)#M\Diagonal(ones(npoin_ang_total))
    
    A = M_inv * LHS
    M = zeros(1)
    M_inv = zeros(1)
    LHS = zeros(1)
    M_sp = zeros(1)
    LHS_sp = zeros(1)
    LHS_el_ang = zeros(1)
    M_rad_el = zeros(1)
    LHS_ang_spat = zeros(1)
    Mass_ang_spat = zeros(1)
    E_rad_el = zeros(1)
    P_rad_el = zeros(1)
    S_rad_el = zeros(1)

    RHS = zeros(TFloat, npoin_ang_total,1)
    ref = zeros(TFloat, npoin_ang_total,1)
        
    for iel=1:nelem
        for i=1:ngl
            for j =1:ngl
                ip = mesh.connijk[iel,i,j]
                x = mesh.x[ip]
                y = mesh.y[ip]
                if (inputs[:adaptive_extra_meshes])

                else

                    for e_ext = 1:extra_mesh.extra_nelem
                        for iθ = 1:extra_mesh.extra_nop[e_ext]+1
                            ip_ext = extra_mesh.extra_connijk[e_ext,iθ]
                            θ = extra_mesh.extra_coords[1,ip_ext]
                            ip_g = (ip-1) * extra_mesh.extra_npoin + ip_ext
                            κip = exp(-((x-3/2)/3)^2)*exp(-y/2)
                            σip = 0.1*κip
                            gip = exp(-((x-3/3)/3)^2)
                            hip = exp(-4*(2-y)/2)
                            sip = exp(-((96/(2*π))*(θ-7*π/5))^2)
                            uip = gip*hip*sip
                            ref[ip_g] = uip
                            propip = cos(θ)*hip*sip*gip*((-2*x)/9 + 2/9) + sin(θ)*gip*sip*2*hip
                            if (ip in mesh.poin_in_bdy_edge)
                                RHS[ip_g] = uip
                                A[ip_g,:] .= 0.0
                                A[ip_g,ip_g] = 1.0
                            else
                                RHS[ip_g] = -(-(user_f!(x,y,θ))*σip + κip*uip +  propip) 
                            end
                        end
                    end
                end
            end
        end
    end
    #A_inv = inv(A)
    @info "built RHS"
    #@info RHS
    @info "solving system"
    solution = solveAx(A, RHS, inputs[:ode_solver])
    @info "done radiation solved"
    @info "dof", npoin_ang_total
    @info "absolute errors, inf, L1, L2", maximum(abs.(solution - ref)), sum(abs.(solution-ref))/npoin_ang_total, sqrt(sum((solution-ref).^2))/npoin_ang_total
    @info "relative errors, inf, L1 ,L2", maximum(abs.(solution - ref))/maximum(abs.(ref)), sum(abs.(solution-ref))/sum(abs.(ref)), sqrt(sum((solution-ref).^2))/sqrt(sum((ref).^2))
    @info "norms", maximum(abs.(solution - ref)), maximum(abs.(ref)), sum(abs.(solution-ref)), sum(abs.(ref)), sqrt(sum((solution-ref).^2)), sqrt(sum((ref).^2))
    A = zeros(1)
    RHS = zeros(1)
    @info "integrating solution and reference in angle"
    int_sol = zeros(TFloat, npoin,1)
    int_ref = zeros(TFloat, npoin,1)
    for iel=1:nelem
        for i=1:ngl
            for j =1:ngl
                ip = mesh.connijk[iel,i,j]
                x = mesh.x[ip]
                y = mesh.y[ip]
                div = 1
                if (i==1 && j == 1) || (i==ngl && j==ngl)
                    if !(ip in  mesh.poin_in_bdy_edge)
                        div = 4
                    end
                elseif (i==1 || j==1 || i==ngl || j==ngl)
                    if !(ip in  mesh.poin_in_bdy_edge)
                        div = 2
                    end 
                end
                for e_ext = 1:extra_mesh.extra_nelem
                    for iθ = 1:extra_mesh.extra_nop[e_ext]+1
                        ip_ext = extra_mesh.extra_connijk[e_ext,iθ]
                        θ = extra_mesh.extra_coords[1,ip_ext]
                        ip_g = (ip-1) * extra_mesh.extra_npoin + ip_ext
                        int_sol[ip] += solution[ip_g]*extra_mesh.ωθ[iθ]/div
                        int_ref[ip] += ref[ip_g]*extra_mesh.ωθ[iθ]/div
                    end
                end
            end
        end
    end
    plot_triangulation(NSD_2D(), mesh, int_ref[:], "ref",  inputs[:output_dir], inputs; iout=1, nvar=1)
    plot_triangulation(NSD_2D(), mesh, int_sol[:], "sol",  inputs[:output_dir], inputs; iout=2, nvar=1)
end

function _rad_element_mass_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, M_rad_el, iel, 
        nop_ang, nelem_ang, Je_ang, ωθ, ψ_ang, QT::Inexact, SD::NSD_2D, AD::ContGal)
    
    for j=1:ngl
        for i=1:ngl
            #ωJac = ω[i]*ω[j]*Je[iel,i,j]
            J = i +(j-1)*(ngl) 
            for e_ext = 1:nelem_ang
                for iθ = 1:nop_ang[e_ext]+1
                    J_ext = iθ
                    ωJac_rad = ωθ[iθ]*Je_ang[e_ext,iθ]
                            
                    for n=1:ngl
                        for m=1:ngl
                            ωJac = ω[m]*ω[n]*Je[iel,m,n]
                            I = m + (n-1)*(ngl)
                            for jθ = 1:nop_ang[e_ext]+1
                                I_ext = jθ
                                M_rad_el[iel,I,J,e_ext,I_ext,J_ext] += ωJac*ωJac_rad*ψ[j,n]*ψ[i,m]*ψ_ang[iθ,jθ]
                            end
                        end
                    end
                end
            end
        end
    end
end

function _rad_element_extinction_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, E_rad_el, iel, κ, nop_ang, nelem_ang, Je_ang, ωθ, ψ_ang, QT::Inexact, SD::NSD_2D, AD::ContGal)

    for j=1:ngl
        for i=1:ngl
            #ωJac = ω[i]*ω[j]*Je[iel,i,j]
            J = i +(j-1)*(ngl)
            for e_ext = 1:nelem_ang
                for iθ = 1:nop_ang[e_ext]+1
                    J_ext = iθ 
                    ωJac_rad = ωθ[iθ]*Je_ang[e_ext,iθ]

                    for n=1:ngl
                        for m=1:ngl
                            ωJac = ω[m]*ω[n]*Je[iel,m,n]
                            I = m + (n-1)*(ngl)
                            for jθ = 1:nop_ang[e_ext]+1
                                I_ext = jθ
                                E_rad_el[iel,I,J,e_ext,I_ext,J_ext] += κ[i,j]*ωJac*ωJac_rad*ψ[j,n]*ψ[i,m]*ψ_ang[iθ,jθ]
                            end
                        end
                    end
                end
            end
        end
    end
end

function _rad_element_internal_propagation_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, dξdx, dηdy, P_rad_el, iel,
        nop_ang, nelem_ang, Je_ang, ωθ, ψ_ang, connijk_ang, coords_ang, QT::Inexact, SD::NSD_2D, AD::ContGal)

    for j=1:ngl
        for i=1:ngl
            #ωJac = ω[i]*ω[j]*Je[iel,i,j]
            J = i +(j-1)*(ngl)
            dξdx_ij = dξdx[iel,i,j]
            dηdy_ij = dηdy[iel,i,j]
            for e_ext = 1:nelem_ang
                for iθ = 1:nop_ang[e_ext]+1
                    J_ext = iθ
                    ωJac_rad = ωθ[iθ]*Je_ang[e_ext,iθ]
                    
                    for n=1:ngl
                        for m=1:ngl
                            ωJac = ω[m]*ω[n]*Je[iel,m,n] 
                            I = m + (n-1)*(ngl)
                            for jθ = 1:nop_ang[e_ext]+1
                                I_ext = jθ
                                i_angle = connijk_ang[e_ext,iθ]

                                θ = coords_ang[1,i_angle]
                                #s = [cos(θ), sin(θ)]
                                P_rad_el[iel,I,J,e_ext,I_ext,J_ext] += (ωJac*ωJac_rad)*(ψ[i,m]*dψ[n,j]*sin(θ)*dηdy_ij + ψ[j,n]*dψ[m,i]*cos(θ)*dξdx_ij)*(ψ_ang[iθ,jθ])
                                #=if (ψ[i,m]*ψ[j,n]*extra_mesh.ψ[iθ,jθ] != 0)
                                    @info P_rad_el[iel,I,J,e_ext,I_ext,J_ext], dψ[k,n]*s[2]*dηdy_ij+dψ[k,m]*s[1]*dξdx_ij, ψ[i,m]*ψ[j,n]*extra_mesh.ψ[iθ,jθ], ωJac*ωJac_rad
                                end=#
                            end
                        end
                    end
                end
            end
        end
    end
end


function _rad_element_scattering_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, S_rad_el, iel, σ, Φ, nop_ang, nelem_ang, Je_ang, ωθ, ψ_ang, connijk_ang, QT::Inexact, SD::NSD_2D, AD::ContGal)

    for j=1:ngl
        for i=1:ngl
            ωJac = ω[i]*ω[j]*Je[iel,i,j]
            J = i +(j-1)*(ngl)
            for e_ext = 1:nelem_ang
                for iθ = 1:nop_ang[e_ext]+1
                    J_ext = iθ 
                    ωJac_rad = ωθ[iθ]*Je_ang[e_ext,iθ]
                    jpθ = connijk_ang[e_ext,iθ]
                    for n=1:ngl
                        for m=1:ngl
                            I = m + (n-1)*(ngl)
                            for jθ = 1:nop_ang[e_ext]+1
                                I_ext = jθ
                                for e_ext_scatter = 1:nelem_ang
                                    for kθ = 1:nop_ang[e_ext]+1
                                        ipθ = connijk_ang[e_ext_scatter,kθ]
                                        ωJac_rad_scatter = ωθ[kθ]*Je_ang[e_ext,kθ]
                                                        
                                        S_rad_el[iel,I,J,e_ext,I_ext,J_ext] += ωJac*ωJac_rad*ωJac_rad_scatter*ψ[j,n]*ψ[i,m]*σ[m,n]*Φ[jpθ,ipθ]*ψ_ang[iθ,jθ]
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function DSS_angular!(M_el_ang, M, nelem, nop_ang, nelem_ang, connijk_ang, ngl, connijk, iel, ::NSD_2D)
    
    for j=1:ngl
        for i=1:ngl
            J = i +(j-1)*(ngl)
            for e_ext = 1:nelem_ang
                for iθ = 1:nop_ang[e_ext]+1
                    J_ext = iθ 
                    JP_ext = connijk_ang[e_ext,iθ]
                    for n=1:ngl
                        for m=1:ngl
                            I  = m + (n-1)*(ngl)
                            for jθ = 1:nop_ang[e_ext]+1
                                I_ext = jθ 
                                IP_ext = connijk_ang[e_ext,jθ]
                                #@info M[iel, I, J]#, M_el_ang[iel, I, J, e_ext, I_ext, J_ext]
                                M[iel, I, J, IP_ext, JP_ext] += M_el_ang[iel, I, J, e_ext, I_ext, J_ext]
                            end
                        end
                    end
                end
            end
        end
    end
end

function DSS_spatial!(M_sp, M_ang_spat, nelem, npoin_ang, ngl, connijk, inputs, ::NSD_2D)

    for iel = 1:nelem
        for j = 1:ngl
            for i=1:ngl
                J = i +(j-1)*(ngl)
                jp = connijk[iel, i, j]
                if (inputs[:adaptive_extra_meshes])

                else
                    for jp_ext = 1:npoin_ang

                        for n=1:ngl
                            for m = 1:ngl
                                I  = m + (n-1)*(ngl)
                                ip = connijk[iel, m, n]

                                for ip_ext = 1:npoin_ang
                                    M_sp[ip, jp, ip_ext, jp_ext] += M_ang_spat[iel, I, J, ip_ext, jp_ext]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function Map_rad_global!(M_sp, M, npoin, npoin_ang, nelem, connijk, ngl, inputs, ::NSD_2D)

    if (inputs[:adaptive_extra_meshes])

    else
        step_ip = 0 
        for ip=1:npoin
            step_jp = 0
            for jp =1:npoin
                npoin_ext_ip = npoin_ang
                npoin_ext_jp = npoin_ang
                for ip_ext = 1:npoin_ext_ip
                    idx_ip = (ip-1)*(npoin_ext_ip) + ip_ext
                    for jp_ext = 1:npoin_ext_jp
                        idx_jp = (jp-1)*(npoin_ext_jp) + jp_ext
                        M[idx_ip, idx_jp] = M_sp[ip, jp, ip_ext, jp_ext]
                    end
                end
                step_jp +=1
            end
            step_ip +=1
        end
    end

end

