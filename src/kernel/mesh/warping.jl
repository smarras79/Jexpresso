function warp_mesh!(mesh,inputs)
  
  am = inputs[:a_mount]
  hm = inputs[:h_mount]
  xc = inputs[:c_mount]
  ztop = maximum(mesh.y)#30000.0
  zsurf = zeros(Float64,mesh.npoin)
  sigma = zeros(Float64,mesh.npoin)
  for i=1:mesh.npoin
    sigma[i] = mesh.y[i]
    
    z = (ztop - zsurf[i])/ztop * sigma[i] + zsurf[i]
    
    mesh.y[i] = z
  end
  if (inputs[:mount_type] == "agnesi")
    am = inputs[:a_mount]
    hm = inputs[:h_mount]
    xc = inputs[:c_mount]
    for ip = 1:mesh.npoin
      x = mesh.x[ip]
      zsurf[ip] = hm/(1+ ((x-xc)/am)^2)
    end
  elseif (inputs[:mount_type] == "schar")
    ac = inputs[:a_mount]
    hc = inputs[:h_mount]
    lambdac = inputs[:lambda_mount]
    for ip = 1:mesh.npoin
      x = mesh.x[ip]
      zsurf[ip] = hc * exp(-(x/ac)^2) * cospi(x/lambdac)^2 
    end
  end

  for ip = 1:mesh.npoin
    sigma[ip] = mesh.y[ip]
    #if (mesh.y[ip] < 10000.0)  
      z = (ztop - zsurf[ip])/ztop * sigma[ip] + zsurf[ip]
      mesh.y[ip] = z
    #=elseif (mesh.y[ip] < 15000.0)
      factor = (15000-mesh.y[ip])/5000.0
      z = (ztop - factor*zsurf[ip])/ztop * sigma[ip] + factor*zsurf[ip]
      mesh.y[ip] = z 
    end=#
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

function warp_mesh_3D!(mesh,inputs)
    if (inputs[:mount_type] == "real topography")
        # find surface heights by reading and interpolating real data onto grid
        fname = inputs[:topo_database]
        fname2 = inputs[:topo_geoid]
        lat_min = inputs[:read_topo_latmin]
        lat_max = inputs[:read_topo_latmax]
        lon_min = inputs[:read_topo_lonmin]
        lon_max = inputs[:read_topo_lonmax]
        zone = inputs[:read_topo_zone]
        xmin = minimum(mesh.x)
        xmax = maximum(mesh.x)
        ymin = minimum(mesh.y)
        ymax = maximum(mesh.y)
        lat, lon, z_topo = extract_region_topography_from_global_data(fname, fname2, lat_max, lon_max, lat_min, lon_min)
        
        x_topo, y_topo = Map_lat_lon_onto_simulation_domain(lat,lon,xmin,xmax,ymin,ymax,zone)
        zsurf = zeros(mesh.npoin)
        
        interpolate_topography_onto_grid!(mesh.x, mesh.y, zsurf, x_topo, y_topo, z_topo)
        ### sigma coordinate topography
        ztop = maximum(mesh.z)
        sigma = zeros(mesh.npoin)
        for ip = 1:mesh.npoin
            sigma[ip] = mesh.z[ip]
            z_new = (ztop - zsurf[ip])/ztop * sigma[ip] + zsurf[ip]
            mesh.z[ip] = z_new
        end

  elseif (inputs[:mount_type] == "agnesi")
    zsurf = zeros(mesh.npoin)
    sigma = zeros(mesh.npoin)
    ztop = maximum(mesh.z)
    am = inputs[:a_mount]
    hm = inputs[:h_mount]
    xc = inputs[:c_mount]
    for ip = 1:mesh.npoin
      x = mesh.x[ip]
      zsurf[ip] = hm/(1+ ((x-xc)/am)^2)
    end
  elseif (inputs[:mount_type] == "schar")
    ac = inputs[:a_mount]
    hc = inputs[:h_mount]
    lambdac = inputs[:lambda_mount]
    for ip = 1:mesh.npoin
      x = mesh.x[ip]
      zsurf[ip] = hc * exp(-(x/ac)^2) * cospi(x/lambdac)^2
    end
  end

  for ip = 1:mesh.npoin
    sigma[ip] = mesh.z[ip]
      z = (ztop - zsurf[ip])/ztop * sigma[ip] + zsurf[ip]
      mesh.z[ip] = z
  end
end


    
