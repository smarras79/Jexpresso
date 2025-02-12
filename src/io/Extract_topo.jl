

function extract_region_topography_from_global_data(fname,fname2, lat_max,lon_max,lat_min,lon_min)

    lat = ncread(fname,"lat")
    lon = ncread(fname,"lon")
    z = ncread(fname, "z")
    geoid = ncread(fname, "z")
    
    z += geoid
    n_lat = size(lat,1)
    n_lon = size(lon,1)
    # find number of points on array
    n_lat_array = 0
    n_lon_array = 0
    for i=1:n_lat
        if (lat[i] <= lat_max && lat[i] >= lat_min)
            n_lat_array +=1 
        end
    end
    for j=1:n_lon
        if (lon[j] <= lon_max && lon[j] >= lon_min)
            n_lon_array +=1
        end
    end
    
    lat_out = zeros(n_lat_array,1)
    lon_out = zeros(n_lon_array,1)
    z_out = zeros(n_lon_array,n_lat_array)

    i_count = 1
    for i=1:n_lat
        if (lat[i] <= lat_max && lat[i] >= lat_min) 
            lat_out[i_count] = lat[i]
            i_count +=1
        end
    end
    j_count=1
    for j=1:n_lon
        if (lon[j] <= lon_max && lon[j] >= lon_min)
            lon_out[j_count] = lon[j]
            j_count +=1
        end
    end
    i_count = 1
    j_count =1
    for i=1:n_lat
        if (lat[i] <= lat_max && lat[i] >= lat_min)

            j_count =1
            for j=1:n_lon
                if (lon[j] <= lon_max && lon[j] >= lon_min)
                    z_out[j_count,i_count] = max(z[j,i],0.0) #ignore bathymetry
                    j_count +=1
                end
            end
            i_count +=1
        end
    end
    
    return lat_out, lon_out, z_out

end

function Map_lat_lon_onto_simulation_domain(lat,lon,xmin,xmax,ymin,ymax,zone)

    ## first project on UTM from WGS84
    n_lat = size(lat,1)
    n_lon = size(lon,1)

    x_utm = zeros(n_lat,n_lon)
    y_utm = zeros(n_lat,n_lon)

    utm_zone = UTMfromLLA(zone, true, wgs84)
    for i=1:n_lat
        for j=1:n_lon
            lla_coord = LLA(lat[i],lon[j])
            utm_coords = utm_zone(lla_coord)
            x_utm[i,j] = utm_coords.x
            y_utm[i,j] = utm_coords.y
        end
    end

    #find the rectangle corners
    utm_x_max = maximum(x_utm)
    utm_x_min = minimum(x_utm)

    utm_y_max = maximum(y_utm)
    utm_y_min = minimum(y_utm)

    ### create a linear mapping between utm rectangle and simulation surface rectangle

    x_factor = (xmax-xmin)/(utm_x_max-utm_x_min)
    y_factor = (ymax-ymin)/(utm_y_max-utm_y_min)
    x = zeros(n_lat,n_lon)
    y = zeros(n_lat,n_lon)

    for i=1:n_lat
        for j=1:n_lon
            x[i,j] = (x_utm[i,j] - utm_x_min)*x_factor + xmin 

            y[i,j] = (y_utm[i,j]-utm_y_min)*y_factor + ymin
        end
    end
    return x,y
end

function interpolate_topography_onto_grid!(x,y,z_surf, x_topo, y_topo, z_topo)

    ymin = minimum(y)
    ymax = maximum(y)
    xmin = minimum(x)
    xmax = maximum(x)
    npoin = size(x,1)
    ni_topo = size(x_topo,1)
    nj_topo = size(x_topo,2)
    for ip=1:npoin
        ### find points for bilinear interpolation
        y11 = ymin
        y12 = ymax
        y21 = ymin
        y22 = ymax
        x11 = xmin
        x21 = xmax
        x12 = xmin
        x22 = xmax
        dist11 =100000000.0
        dist12 =100000000.0
        dist21 =100000000.0
        dist22 =100000000.0

        z11 = 0.0
        z12 = 0.0
        z21 = 0.0
        z22 = 0.0

        for i = 1:ni_topo
            for j=1:nj_topo
                dist = sqrt((x[ip]-x_topo[i,j])^2+(y[ip]-y_topo[i,j])^2)
                if (x[ip] >= x_topo[i,j] && y[ip] >= y_topo[i,j] && dist <=dist11)
                    x11 = x_topo[i,j]
                    y11 = y_topo[i,j]
                    z11 = z_topo[j,i]
                    dist11 = dist
                end
                if (x[ip] >= x_topo[i,j] && y[ip] <= y_topo[i,j] && dist <=dist12)
                    x12 = x_topo[i,j]
                    y12 = y_topo[i,j]
                    z12 = z_topo[j,i]
                    dist12 = dist
                end
                if (x[ip] <= x_topo[i,j] && y[ip] >= y_topo[i,j] && dist <=dist21)
                    x21 = x_topo[i,j]
                    y21 = y_topo[i,j]
                    z21 = z_topo[j,i]
                    dist21 = dist
                end
                if (x[ip] <= x_topo[i,j] && y[ip] <= y_topo[i,j] && dist <=dist22)
                    x22 = x_topo[i,j]
                    y22 = y_topo[i,j]
                    z22 = z_topo[j,i]
                    dist22 = dist
                end
            end
        end
        y1 = max(y11,y21)
        y2 = min(y22,y12)
        x1 = max(x11,x12)
        x2 = min(x22,x21)

        #@info "x",x[ip],x11,x12,x21,x22
        #@info "y",y[ip],y11,y12,y21,y22
        ## perform bilinear interpolation
        if (y2 > y1 && x2 > x1) #&& (y[ip] > ymin && y[ip] > ymax && x[ip] > xmin && x[ip] < xmax)
            t1 = (y2 - y[ip])/(y2-y1)*((x2-x[ip])/(x2-x1)*z11+ (x[ip]-x1)/(x2-x1)*z21)
            t2 = (y[ip]-y1)/(y2-y1)*((x2-x[ip])/(x2-x1)*z12+ (x[ip]-x1)/(x2-x1)*z22)
            z_surf[ip] = t1 + t2
        else
            z_surf[ip] = 0.0
        end
        #elseif(y2 > y1)
        #    z_surf[ip] = z11 + (y[ip] - y1)*(z12-z11)/(y2-y1)
        #elseif(x2 > x1)
        #    z_surf[ip] = z11 + (x[ip] - x1)*(z21-z11)/(x2-x1)
        #else
        #    z_surf[ip] = z11 
        #end
    end
end
