function _rad_element_mass_matrix_3Dby2D!(u, neqs, ngl, dψ, ψ, ω, Je, M_rad_el, iel, 
        extra_mesh::St_extra_mesh,::CL, QT::Inexact, SD::NSD_3D, AD::ContGal)
    
    for k=1:ngl
        for j=1:ngl
            for i=1:ngl
                ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                J = i +(j-1)*(ngl) + k*(ngl)*(ngl) 
                for e_ext = 1:extra_mesh.nelem
                    for iϕ = 1:extra_mesh.nop[e_ext]+1
                        for iθ = 1:extra_mesh.nop[e_ext]+1
                            J_ext = iθ + iϕ*(extra_mesh.nop[e_ext]+1)
                            ωJac_rad = extra_mesh.ωθ[iθ]*extra_mesh.ωϕ[iϕ]*extra_mesh.extra_metrics.Je[e_ext,iθ,iϕ]
                            
                            for o=1:ngl
                                for n=1:ngl
                                    for m=1:ngl
                                        I = m + (n-1)*(ngl) + (o-1)*(ngl)*(ngl)
                                        for jϕ = 1:extra_mesh.nop[e_ext]+1
                                            for jθ = 1:extra_mesh.nop[e_ext]+1
                                                I_ext = jθ + jϕ*(extra_mesh.nop[e_ext]+1)
                                                M_rad_el[iel,I,J,e_ext,I_ext,J_ext] += ωJac*ωJac_rad*ψ[o,k]*ψ[j,n]*ψ[i,m]*extra_mesh.basis.ψ[iϕ,jϕ]*extra_mesh.basis.ψ[iθ,jθ]
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
    end
end

function _rad_element_extinction_matrix_3Dby2D!(u, neqs, ngl, dψ, ψ, ω, Je, E_rad_el, iel, κ
        extra_mesh::St_extra_mesh,::CL, QT::Inexact, SD::NSD_3D, AD::ContGal)

    for k=1:ngl
        for j=1:ngl
            for i=1:ngl
                ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                J = i +(j-1)*(ngl) + k*(ngl)*(ngl)
                for e_ext = 1:extra_mesh.nelem
                    for iϕ = 1:extra_mesh.nop[e_ext]+1
                        for iθ = 1:extra_mesh.nop[e_ext]+1
                            J_ext = iθ + iϕ*(extra_mesh.nop[e_ext]+1)
                            ωJac_rad = extra_mesh.ωθ[iθ]*extra_mesh.ωϕ[iϕ]*extra_mesh.extra_metrics.Je[e_ext,iθ,iϕ]

                            for o=1:ngl
                                for n=1:ngl
                                    for m=1:ngl
                                        I = m + (n-1)*(ngl) + (o-1)*(ngl)*(ngl)
                                        for jϕ = 1:extra_mesh.nop[e_ext]+1
                                            for jθ = 1:extra_mesh.nop[e_ext]+1
                                                I_ext = jθ + jϕ*(extra_mesh.nop[e_ext]+1)
                                                E_rad_el[iel,I,J,e_ext,I_ext,J_ext] += κ[i,j,k]*ωJac*ωJac_rad*ψ[o,k]*ψ[j,n]*ψ[i,m]*extra_mesh.basis.ψ[iϕ,jϕ]*extra_mesh.basis.ψ[iθ,jθ]
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
    end
end

function _rad_element_internal_propagation_matrix_3Dby2D!(u, neqs, ngl, dψ, ψ, ω, Je, P_rad_el, iel,
        extra_mesh::St_extra_mesh,::CL, QT::Inexact, SD::NSD_3D, AD::ContGal)

    for k=1:ngl
        for j=1:ngl
            for i=1:ngl
                ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                J = i +(j-1)*(ngl) + k*(ngl)*(ngl)
                for e_ext = 1:extra_mesh.nelem
                    for iϕ = 1:extra_mesh.nop[e_ext]+1
                        for iθ = 1:extra_mesh.nop[e_ext]+1
                            J_ext = iθ + iϕ*(extra_mesh.nop[e_ext]+1)
                            ωJac_rad = extra_mesh.ωθ[iθ]*extra_mesh.ωϕ[iϕ]*extra_mesh.extra_metrics.Je[e_ext,iθ,iϕ]

                            for o=1:ngl
                                for n=1:ngl
                                    for m=1:ngl
                                        I = m + (n-1)*(ngl) + (o-1)*(ngl)*(ngl)
                                        for jϕ = 1:extra_mesh.nop[e_ext]+1
                                            for jθ = 1:extra_mesh.nop[e_ext]+1
                                                I_ext = jθ + jϕ*(extra_mesh.nop[e_ext]+1)
                                                i_angle = extra_mesh.extra_connijk[e_ext,jθ,jϕ]

                                                θ = extra_mesh.extra_coords[1,i_angle]
                                                ϕ = extra_mesh.extra_coords[2,i_angle]
                                                s = [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]
                                                P_rad_el[iel,I,J,e_ext,I_ext,J_ext] += (ωJac*ωJac_rad)*(ψ[k,o]*dψ[k,o]*s[3]+ψ[n,j]*dψ[j,n]*s[2]+
                                                                                                       ψ[i,m]*dψ[i,m]*s[1])*(extra_mesh.basis.ψ[iϕ,jϕ]*extra_mesh.basis.ψ[iθ,jθ])
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
    end
end
   #=                                                                    
                            #=dξdx_ij = dξdx[iel,i,j,k]
                            dξdy_ij = dξdy[iel,i,j,k]
                            dξdz_ij = dξdz[iel,i,j,k]

                            dηdx_ij = dηdx[iel,i,j,k]
                            dηdy_ij = dηdy[iel,i,j,k]
                            dηdz_ij = dηdz[iel,i,j,k]

                            dζdx_ij = dζdx[iel,i,j,k]
                            dζdy_ij = dζdy[iel,i,j,k]
                            dζdz_ij = dζdz[iel,i,j,k]

                            dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij + dFdζ*dζdx_ij
                            dGdx = dGdξ*dξdx_ij + dGdη*dηdx_ij + dGdζ*dζdx_ij
                            dHdx = dHdξ*dξdx_ij + dHdη*dηdx_ij + dHdζ*dζdx_ij

                            dFdy = dFdξ*dξdy_ij + dFdη*dηdy_ij + dFdζ*dζdy_ij
                            dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij + dGdζ*dζdy_ij
                            dHdy = dHdξ*dξdy_ij + dHdη*dηdy_ij + dHdζ*dζdy_ij

                            dFdz = dFdξ*dξdz_ij + dFdη*dηdz_ij + dFdζ*dζdz_ij
                            dGdz = dGdξ*dξdz_ij + dGdη*dηdz_ij + dGdζ*dζdz_ij
                            dHdz = dHdξ*dξdz_ij + dHdη*dηdz_ij + dHdζ*dζdz_ij
                            =#
                            s = [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)] 
                            auxi = ωJac_rad*((s[1]*dFdx + s[2]*dGdy + s[3]*dHdz) - S[i,j,k,e_ext,iθ,iϕ])
                            rhs_extra_el[e_ext,iθ,iϕ] -= auxi
                        end
                    end
                end
                DSS_extra!(extra_mesh, rhs_extra_el, RHS_extra,iel,i,j,k)
                for ip_extra = 1:extra_mesh.npoin
                    auxi = ωJac*RHS_extra[i,j,k,ip_extra]
                    rhs_el[iel,i,j,k,ip] -= auxi
                end
            end
        end
    end
end=#
function _rad_element_scattering_matrix_3Dby2D!(u, neqs, ngl, dψ, ψ, ω, Je, S_rad_el, iel,σ,Φ
        extra_mesh::St_extra_mesh,::CL, QT::Inexact, SD::NSD_3D, AD::ContGal)

    for k=1:ngl
        for j=1:ngl
            for i=1:ngl
                ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                J = i +(j-1)*(ngl) + k*(ngl)*(ngl)
                for e_ext = 1:extra_mesh.nelem
                    for iϕ = 1:extra_mesh.nop[e_ext]+1
                        for iθ = 1:extra_mesh.nop[e_ext]+1
                            J_ext = iθ + iϕ*(extra_mesh.nop[e_ext]+1)
                            ωJac_rad = extra_mesh.ωθ[iθ]*extra_mesh.ωϕ[iϕ]*extra_mesh.extra_metrics.Je[e_ext,iθ,iϕ]

                            for o=1:ngl
                                for n=1:ngl
                                    for m=1:ngl
                                        I = m + (n-1)*(ngl) + (o-1)*(ngl)*(ngl)
                                        for jϕ = 1:extra_mesh.nop[e_ext]+1
                                            for jθ = 1:extra_mesh.nop[e_ext]+1
                                                I_ext = jθ + jϕ*(extra_mesh.nop[e_ext]+1)
                                                for e_ext_scatter = 1:extra_mesh.nelem
                                                    for kϕ = 1:extra_mesh.nop[e_ext_scatter]+1
                                                        for kθ = 1:extra_mesh.nop[e_ext_scatter]+1
                                                            ωJac_rad_scatter = extra_mesh.ωθ[kθ]*extra_mesh.ωϕ[kϕ]*extra_mesh.extra_metrics.Je[e_ext,kθ,kϕ]
                                                        
                                                            S_rad_el[iel,I,J,e_ext,I_ext,J_ext] += ωJac*ωJac_rad*ωJac_rad_scatter*ψ[o,k]*ψ[j,n]*ψ[i,m]*
                                                            extra_mesh.basis.ψ[iϕ,jϕ]*extra_mesh.basis.ψ[iθ,jθ]*σ[m,n,o]*Φ[jθ,jϕ,kθ,kϕ]
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
                end
            end
        end
    end
end

function DSS_spatial_angular(M_el_ang, M, nelem, extra_mesh, ngl, connijk)
    
    for e=1:nelem
        for k=1:ngl
            for j=1:ngl
                for i=1:ngl
                    J = i +(j-1)*(ngl) + k*(ngl)*(ngl)
                    JP = connijk[e,i,j,k]
                    for e_ext = 1:extra_mesh.nelem
                        for iϕ = 1:extra_mesh.nop[e_ext]+1
                            for iθ = 1:extra_mesh.nop[e_ext]+1
                                J_ext = iθ + iϕ*(extra_mesh.nop[e_ext]+1)
                                JP_ext = extra_mesh.extra_connijk[e_ext,iθ,iϕ]
                                for o=1:ngl
                                    for n=1:ngl
                                        for m=1:ngl
                                            I  = m + (n-1)*(ngl) + (o-1)*(ngl)*(ngl)
                                            IP = connijk[e,m,n,o]
                                            for jϕ = 1:extra_mesh.nop[e_ext]+1
                                                for jθ = 1:extra_mesh.nop[e_ext]+1
                                                    I_ext = jθ + jϕ*(extra_mesh.nop[e_ext]+1)
                                                    IP_ext = extra_mesh.extra_connijk[e_ext,jθ,jϕ]
                                                    M[IP,JP].M_ang[IP_ext,JP_ext] = M[IP,JP].M_ang[IP_ext,JP_ext] + M_el_ang[e,I,J,e_ext,I_ext,J_ext]
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
        end
    end
end

function Map_rad_global(M_ang_spat, M, npoin,extra_meshes)
    step_ip = 0.0 
    for ip=1:npoin
        step_jp = 0.0 
        for jp =1:npoin
            npoin_ext_ip = extra_meshes[ip].npoin
            npoin_ext_jp = extra_meshes[ip].npoin
            for ip_ext = 1:npoin_ext_ip
                idx_ip = ip*(step_ip+1) + ip_ext-1
                for jp_ext = 1:npoin_ext_jp
                    idx_jp = jp*(step_jp+1) + jp_ext-1
                    M[idx_ip,idx_jp] = M_ang_spat[ip,jp].M_ang[ip_ext,jp_ext]
                     
                    

                end
            end
            step_jp +=1
        end
        step_ip +=1
    end

end

