function user_source(q::Array, mesh::St_mesh, T)

    S = zeros(T, mesh.npoin)

    #
    # S(x,y) .= sinpi.(c*(x .- xc)).*sinpi(c*(y .- yc))    
    #
    c = 2.0
    xc, yc = (maximum(mesh.x) + minimum(mesh.x))/2, (maximum(mesh.y) + minimum(mesh.y))/2
    for iel_g = 1:mesh.nelem
        for i=1:mesh.ngl
            for j=1:mesh.ngl
                
                ip = mesh.connijk[i,j,iel_g];
                x, y = mesh.x[ip], mesh.y[ip];
                
                S[ip] = -2.0*(c*π)^2*sinpi(c*(x - xc))*sinpi(c*(y - yc))
                
            end
        end
    end
    
    return  S
    
end
