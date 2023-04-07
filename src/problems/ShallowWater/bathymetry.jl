function bathymetry(x,y)

    return -2.0

end

#=function bathymetry(x)

    return -2.0

end=#

#1D parabolic bump Swashes case 1
function bathymetry(x)

    if (x >= 8 && x <= 12)
        Hb = 0.2 - 0.05 * (x - 10)^2
    else
        Hb = 0.0
    end 
    return Hb
end

function compute_variable_bathymetry(SD::NSD_1D, QT::Inexact, q, mesh, metrics, basis, ω, M, T)
    H = zeros(mesh.ngl,mesh.nelem)
    dHdx_el = zeros(mesh.ngl,mesh.nelem)
    dzdx_el = zeros(mesh.ngl,mesh.nelem)
    for e=1:mesh.nelem
        for i=1:mesh.ngl
            ip = conn[i,e]
            H[i,e] = q[ip,1]
        end
        for i=1:mesh.ngl
            dHdξ = 0.0
            for k=1:mesh.ngl
                dHdξ  = dHdξ + basis.dψ[k,i]*H[k,iel]
            end
            dHdx_el[i,e] += dHdξ
        end
    end
    for e=1:mesh.nelem
        for i=1:mesh.ngl
            ip = conn[i,e]
            dzdx_el[i,e] = (ω[i]*mesh.Δx[iel]/2)*(((q[ip,2]^2)/(9.81*H[i,e]^3))*dHdx_el[i,e] )#Add friction term here
        end
    end
    dzdx = DSS_rhs(SD, dzdx_el, mesh.connijk, mesh.nelem, mesh.npoin, 1, mesh.nop, T)
    divive_by_mass_matrix!(dzdx, M, QT,1)

    return dzdx
end




