function warp_mesh!(mesh,inputs)
  
  am = 1000.0
  hm = 1.0
  xc = 0.0
  ztop = 30000.0
  zsurf = zeros(Float64,mesh.npoin)
  sigma = zeros(Float64,mesh.npoin)
  for i=1:mesh.npoin
    sigma[i] = mesh.y[i]
    
    z = (ztop - zsurf[i])/ztop * sigma[i] + zsurf[i]
    
    mesh.y[i] = z
  end

  for ip = 1:mesh.npoin
    x = mesh.x[ip]
    zsurf[ip] = hm/(1+ ((x-xc)/am)^2)
  end

  for ip = 1:mesh.npoin
    sigma[ip] = mesh.y[ip] 
    z = (ztop - zsurf[ip])/ztop * sigma[ip] + zsurf[ip]
    mesh.y[ip] = z
  end 
  
  #=for iedge = 1:size(mesh.bdy_edge_in_elem,1)
        iel = mesh.bdy_edge_in_elem[iedge]
        comp = mesh.bdy_edge_comp[iedge]
        for k=1:mesh.ngl
            #if (mesh.bdy_edge_type[iedge] == "free_slip")
                tag = mesh.bdy_edge_type[iedge]
                ip = mesh.poin_in_bdy_edge[iedge,k]
                @info "prewarp", mesh.y[ip]
                if (mesh.y[ip] < 10.0)
                  mesh.y[ip] = (hm*(am^2))/((mesh.x[ip]-xc)^2 + am^2)
                end
                @info "postwarp", mesh.y[ip]
            #end
        end
   end=#

end


