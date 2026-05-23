
function restructure4periodicity_1D!(mesh,
                                     coords,
                                     xmax,xmin,ymax,ymin,zmax,zmin,
                                     poin_in_bdy_face,poin_in_bdy_edge,ngl,ngr,nelem,npoin,nsd,bdy_edge_type,
                                     bdy_face_type,bdy_face_in_elem,bdy_edge_in_elem,
                                     connijk,connijk_lag,npoin_linear,nelem_semi_inf,
                                     inputs,backend)
    
    #determine boundary vectors
    x_spare = zeros(npoin,1)
    y_spare = zeros(npoin,1) 
    z_spare = zeros(npoin,1)
    connijk_spare = zeros(nelem,ngl,ngl,ngl)
    poin_in_bdy_face_spare = zeros(size(poin_in_bdy_face,1),ngl,ngl)
    if (nsd == 2)
	if (inputs[:lperiodic_laguerre] && "Laguerre" in bdy_edge_type)
            xmin = minimum(coords[:,1])
            xmax = maximum(coords[:,1])
            e1 = 0
            e = 1
            i1 = 0
            i2 = 0
            while (e1 == 0)
                ip_temp = connijk_lag[e,1,1]
                ip_temp1 = connijk_lag[e,ngl,1]
                if (AlmostEqual(coords[ip_temp,1],xmin))
                    e1 = e
                    i1 = 1
                end
                if (AlmostEqual(coords[ip_temp1,1],xmin))
                    e1 = e
                    i1 = ngl
                end
                e += 1

            end
            e2 = 0
            e = 1
            while (e2 == 0)
                ip_temp = connijk_lag[e,1,1]
                ip_temp1 = connijk_lag[e,ngl,1]
                if (AlmostEqual(coords[ip_temp,1],xmax))
                    e2 = e
                    i2 = 1
                end
                if (AlmostEqual(coords[ip_temp1,1],xmax))
                    e2 = e
                    i2 = ngl
                end
                e += 1

            end
            for j = 1:ngr
                ip1 = connijk_lag[e1,i1,j]
                ip2 = connijk_lag[e2,i2,j]
                if (coords[ip1,1] < coords[ip2,1])
                    ip_dest = ip1
                    ip_kill = ip2
                elseif(coords[ip1,1] > coords[ip2,1])
                    ip_dest = ip2
                    ip_kill = ip1
                end
                connijk_lag[e1,i1,j] = ip_dest
                connijk_lag[e2,i2,j] = ip_dest
	        if (j > 1)
                    for ip = ip_kill:npoin-1
                        coords[ip,1] = coords[ip+1,1]
                        coords[ip,2] = coords[ip+1,2]
                    end
                    for e=1:nelem_semi_inf
                        for i=1:ngl
                            for j=1:ngr
                                ip = connijk_lag[e,i,j]
                                if (ip > ip_kill)
                                    connijk_lag[e,i,j] = ip-1
                                end
                            end
                        end
                    end
                else
                    #  coords[ip_dest,1] = xmin
                end
            end
            npoin -= (ngr-1)
        end  
        if ("periodicx" in bdy_edge_type)
            finder = false
            iedge_bdy = 1
            while (finder == false)
                if (bdy_edge_type[iedge_bdy] == "periodicx")
                    ip = poin_in_bdy_edge[iedge_bdy,1]
                    ip1 = poin_in_bdy_edge[iedge_bdy,2]
                    per1 = [coords[ip,1] - coords[ip1,1], coords[ip,2] - coords[ip1,2]]
                    finder = true
                else
                    iedge_bdy +=1
                end
            end
        else
            per1 = [0.0, 1.0]
        end
        if ("periodicz" in bdy_edge_type)
            finder = false
            iedge_bdy = 1
            while (finder == false)
                if (bdy_edge_type[iedge_bdy] == "periodicz")
                    ip = poin_in_bdy_edge[iedge_bdy,1]
                    ip1 = poin_in_bdy_edge[iedge_bdy,2]
                    per2 = [coords[ip,1] - coords[ip1,1], coords[ip,2] - coords[ip1,2]]
                    finder = true
                else
                    iedge_bdy +=1
                end
            end 
        else
            per2 = [1.0, 0.0]
        end     
        xx = zeros(TFloat, size(coords[:,1],1),1)#KernelAbstractions.zeros(CPU(), TFloat, size(x,1),1)
        yy = zeros(TFloat, size(coords[:,2],1),1)#KernelAbstractions.zeros(CPU(), TFloat, size(y,1),1)
        poin_bdy= zeros(TInt, size(poin_in_bdy_edge))#KernelAbstractions.zeros(CPU(), TInt,size(poin_in_bdy_edge))
        xx .= coords[:,1]
        yy .= coords[:,2]
        poin_bdy .=poin_in_bdy_edge
        interval = [2,3,4]
        for iedge_bdy =1:size(bdy_edge_type,1)
            if ("periodicz" in bdy_edge_type && "periodicx" in bdy_edge_type)
                iel = bdy_edge_in_elem[iedge_bdy]
                for k=1:ngl
                    ip = poin_in_bdy_edge[iedge_bdy,k]
                    if (ip in interval)
                        ip_kill=ip
                        ip_dest=1
                        for e=1:nelem
                            for ii=1:ngl
                                for jj=1:ngl
                                    if (connijk[iel,ii,jj] == ip_kill)
                                        connijk[iel,ii,jj] = ip_dest
                                    end
                                end
                            end
                        end
                        for i=ip_kill:npoin-1
                            coords[i,1] = coords[i+1,1]
                            coords[i,2] = coords[i+1,2]
                        end
                        npoin = npoin-1
                        for iedge =1:size(poin_in_bdy_edge,1)
                            for kk=1:ngl
                                val = poin_in_bdy_edge[iedge,kk]
                                if (val > ip_kill)
                                    poin_in_bdy_edge[iedge,kk] = val - 1
                                end
                                if (val == ip_kill)
                                    poin_in_bdy_edge[iedge,kk] = ip_dest
                                end
                            end
                        end

                        for e=1:nelem
                            for ii=1:ngl
                                for jj=1:ngl
                                    ipp = connijk[e,ii,jj]
                                    if (ipp > ip_kill)
                                        connijk[e,ii,jj] -= 1
                                    end
                                end
                            end
                        end
                        if (inputs[:lperiodic_laguerre])
                            for e=1:nelem_semi_inf
                                for ii=1:ngl
                                    for jj=1:ngr
                                        ipp = connijk_lag[e,ii,jj]
                                        if (ipp > ip_kill)
                                            connijk_lag[e,ii,jj] -= 1
                                        end
                                    end
                                end
                            end
                        end
                        for i=1:size(interval,1)
                            if (interval[i] > ip_kill)
                                interval[i] -= 1
                            end
                            if (interval[i] == ip_kill)
                                interval[i] = 0
                            end 
                        end
                    end   
                end
            end
        end
        # New periodicity interface
        for iedge_bdy =1:size(bdy_edge_type,1)
            if (bdy_edge_type[iedge_bdy] == "periodicx" || bdy_edge_type[iedge_bdy] == "periodicz")
                iel = bdy_edge_in_elem[iedge_bdy]
                for k =1:ngl
                    ip = poin_bdy[iedge_bdy, k]
                    ip_true = poin_in_bdy_edge[iedge_bdy,k]
                    x1 = coords[ip,1]
                    y1 = coords[ip,2]
                    m=1
                    l=1
                    for ii=1:ngl
                        for jj=1:ngl
                            if (connijk[iel,ii,jj] == ip_true)
                                l=ii
                                m=jj
                            end
                        end
                    end
                    if (k < ngl)
                        ip1 = poin_bdy[iedge_bdy,k+1]
                    else
                        ip1 = poin_bdy[iedge_bdy,k-1]
                    end
                    x3 = coords[ip1,1]
                    y3 = coords[ip1,2]
                    vec_bdy = [x1-x3,y1-y3]
                    #@info vec_bdy, per1, per2, determine_colinearity(vec_bdy,per2), determine_colinearity(vec_bdy,per1),x1, x3, y1, y3
                    if (determine_colinearity(vec_bdy,per1))
                        per = per2
                    elseif (determine_colinearity(vec_bdy,per2))
                        per = per1
                    else
                        error("periodicity requested but boundaries cannot match any vectors")
                    end
                    # find corresponding periodic edge point
                    for iedge_per = iedge_bdy+1:size(bdy_edge_type,1)
                        if (bdy_edge_type[iedge_per] == bdy_edge_type[iedge_bdy])
                            iel_per = bdy_edge_in_elem[iedge_per]
                            for k_per=1:ngl
                                ip_per = poin_bdy[iedge_per,k_per]
                                ip_true1 = poin_in_bdy_edge[iedge_per,k_per]
                                x2 = xx[ip_per]
                                y2 = yy[ip_per]
                                if (k_per < ngl)
                                    ip_per1 = poin_bdy[iedge_per,k_per+1]
                                else
                                    ip_per1 = poin_bdy[iedge_per,k_per-1]
                                end
                                x4 = xx[ip_per1]
                                y4 = yy[ip_per1]
                                vec_per = [x2-x4, y2-y4]

                                vec = [x1 - x2, y1 - y2]
                                #check colinearity with periodicity vectors
                                #@info "looking for match", vec,per,vec_bdy,vec_per,determine_colinearity(vec,per), determine_colinearity(vec_per,vec_bdy)

                                if (determine_colinearity(vec,per) && determine_colinearity(vec_per,vec_bdy))
                                    m1=1
                                    l1=1
                                    for ii=1:ngl
                                        for jj=1:ngl
                                            if (connijk[iel_per,ii,jj] == ip_true1)
                                                l1=ii
                                                m1=jj
                                            end
                                        end
                                    end

                                    if (coords[ip_true,1] < coords[ip_true1,1])
                                        ip_dest = ip_true
                                        ip_kill = ip_true1
                                    elseif (coords[ip_true,1] >= coords[ip_true1,1])
                                        ip_dest = ip_true1
                                        ip_kill = ip_true
                                    end
                                    if (coords[ip_true,1]*abs(coords[ip_true,2]) < coords[ip_true1,1]*abs(coords[ip_true1,2]) || coords[ip_true,2]*abs(coords[ip_true,1]) < coords[ip_true1,2]*abs(coords[ip_true1,1]) )
                                        ip_dest = ip_true
                                        ip_kill = ip_true1
                                    else
                                        ip_dest = ip_true1
                                        ip_kill = ip_true
                                    end 
                                    connijk[iel_per,l1,m1] = ip_dest
                                    connijk[iel,l,m] = ip_dest
                                    poin_in_bdy_edge[iedge_per,k_per] = ip_dest
                                    poin_in_bdy_edge[iedge_bdy,k] = ip_dest
                                    ip_true = ip_dest
                                    if !(ip_kill in connijk)

                                        for i=ip_kill:npoin-1
                                            coords[i,1] = coords[i+1,1]
                                            coords[i,2] = coords[i+1,2]
                                        end
                                        npoin = npoin-1
                                        for iedge =1:size(poin_in_bdy_edge,1)
                                            for kk=1:ngl
                                                val = poin_in_bdy_edge[iedge,kk]
                                                if (val > ip_kill)
                                                    poin_in_bdy_edge[iedge,kk] = val - 1
                                                end
                                                if (val == ip_kill)
                                                    poin_in_bdy_edge[iedge,kk] = ip_dest
                                                end
                                            end
                                        end

                                        for e=1:nelem
                                            for ii=1:ngl
                                                for jj=1:ngl
                                                    ipp = connijk[e,ii,jj]
                                                    if (ipp > ip_kill)
                                                        connijk[e,ii,jj] -= 1
                                                    end
                                                    if (ipp == ip_kill)
                                                        connijk[e,ii,jj] = ip_dest
                                                    end
                                                end
                                            end
                                        end
                                        if (inputs[:lperiodic_laguerre])
                                            for e=1:nelem_semi_inf
                                                for ii=1:ngl
                                                    for jj=1:ngr
                                                        ipp = connijk_lag[e,ii,jj]
                                                        if (ipp > ip_kill)
                                                            connijk_lag[e,ii,jj] -= 1
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
        mesh.npoin = npoin
    elseif (nsd == 3)
        if ("periodicx" in bdy_face_type)
            finder = false
            iface_bdy = 1
            while (finder == false)
                if (bdy_face_type[iface_bdy] == "periodicx")
                    ip = poin_in_bdy_face[iface_bdy,1,1]
                    ip1 = poin_in_bdy_face[iface_bdy,1,2]
                    ip2 = poin_in_bdy_face[iface_bdy,2,1]
                    t1 = [coords[ip,1] - coords[ip1,1], coords[ip,2] - coords[ip1,2], coords[ip,3] - coords[ip1,3]]
                    t2 = [coords[ip,1] - coords[ip2,1], coords[ip,2] - coords[ip2,2], coords[ip,3] - coords[ip2,3]]
                    s1 = t1[2]*t2[3] - t1[3]*t2[2]
                    s2 = t1[3]*t2[1] - t1[1]*t2[3]
                    s3 = t1[1]*t2[2] - t1[2]*t2[1]
                    mag = sqrt(s1^2 + s2^2 + s3^2)
                    nor1 = [s1/mag, s2/mag, s3/mag]
                    finder = true
                else
                    iface_bdy +=1
                end
            end
        else
            nor1 = [1.0, 0.0, 0.0]
        end
        if ("periodicz" in bdy_face_type)
            finder = false
            iface_bdy = 1
            while (finder == false)
                if (bdy_face_type[iface_bdy] == "periodicz")
                    ip = poin_in_bdy_face[iface_bdy,1,1]
                    ip1 = poin_in_bdy_face[iface_bdy,1,2]
                    ip2 = poin_in_bdy_face[iface_bdy,2,1]
                    t1 = [coords[ip,1] - coords[ip1,1],coords[ip,2] - coords[ip1,2], coords[ip,3] - coords[ip1,3]]
                    t2 = [coords[ip,1] - coords[ip2,1],coords[ip,2] - coords[ip2,2], coords[ip,3] - coords[ip2,3]]
                    s1 = t1[2]*t2[3] - t1[3]*t2[2]
                    s2 = t1[3]*t2[1] - t1[1]*t2[3]
                    s3 = t1[1]*t2[2] - t1[2]*t2[1]
                    mag = sqrt(s1^2 + s2^2 + s3^2)
                    nor2 = [s1/mag, s2/mag, s3/mag] 
                    finder = true
                else
                    iface_bdy +=1
                end
            end
        else
            nor2 = [0.0, 1.0, 0.0]
        end
        if ("periodicy" in bdy_face_type)
            finder = false
            iface_bdy = 1
            while (finder == false)
                if (bdy_face_type[iface_bdy] == "periodicy")
                    ip = poin_in_bdy_face[iface_bdy,1,1]
                    ip1 = poin_in_bdy_face[iface_bdy,1,2]
                    ip2 = poin_in_bdy_face[iface_bdy,2,1]
                    t1 = [coords[ip,1] - coords[ip1,1], coords[ip,2] - coords[ip1,2], coords[ip,3] - coords[ip1,3]]
                    t2 = [coords[ip,1] - coords[ip2,1], coords[ip,2] - coords[ip2,2], coords[ip,3] - coords[ip2,3]]
                    s1 = t1[2]*t2[3] - t1[3]*t2[2]
                    s2 = t1[3]*t2[1] - t1[1]*t2[3]
                    s3 = t1[1]*t2[2] - t1[2]*t2[1]
                    mag = sqrt(s1^2 + s2^2 + s3^2)
                    nor3 = [s1/mag, s2/mag, s3/mag]
                    finder = true
                else
                    iface_bdy +=1
                end
            end
        else
            nor3 = [0.0, 0.0, 1.0]
        end

        ### triple periodicity for corners 
        interval = [2,3,4,5,6,7,8]
        if ("periodicx" in bdy_face_type && "periodicz" in bdy_face_type && "periodicy" in bdy_face_type)    
            for iface_bdy =1:size(bdy_face_type,1)
                iel = bdy_face_in_elem[iface_bdy]
                for k=1:ngl
                    for l=1:ngl
                        ip = poin_in_bdy_face[iface_bdy,k,l]
                        if (ip in interval)
                            ip_kill=ip
                            ip_dest=1
                            for e=1:nelem
                                for ii=1:ngl
                                    for jj=1:ngl
                                        for kk=1:ngl
                                            if (connijk[iel,ii,jj,kk] == ip_kill)
                                                connijk[iel,ii,jj,kk] = ip_dest
                                            end
                                        end
                                    end
                                end
                            end
                            for i=ip_kill:npoin-1
                                coords[i,1] = coords[i+1,1]
                                coords[i,2] = coords[i+1,2]
                                coords[i,3] = coords[i+1,3]
                            end
                            npoin = npoin-1
                            for iface =1:size(poin_in_bdy_face,1)
                                for kk=1:ngl
                                    for ll=1:ngl
                                        val = poin_in_bdy_face[iface,kk,ll]
                                        if (val > ip_kill)
                                            poin_in_bdy_face[iface,kk,ll] = val - 1
                                        end
                                        if (val == ip_kill)
                                            poin_in_bdy_face[iface,kk,ll] = ip_dest
                                        end
                                    end
                                end
                            end

                            for e=1:nelem
                                for ii=1:ngl
                                    for jj=1:ngl
                                        for kk=1:ngl
                                            ipp = connijk[e,ii,jj,kk]
                                            if (ipp > ip_kill)
                                                connijk[e,ii,jj,kk] -= 1
                                            end
                                        end
                                    end
                                end
                            end
                            for i=1:size(interval,1)
                                if (interval[i] > ip_kill)
                                    interval[i] -= 1
                                elseif (interval[i] == ip_kill)
                                    interval[i] = 0
                                end
                            end
                        end
                    end
                end
            end

            ### done with triple periodicity
        end
### double periodicity
### check if we have more than one periodic set of boundaries
nperiodic = 0
double1 = [0, 0, 0]
double2 = [0, 0, 0]
double3 = [0, 0, 0]
if ("periodicx" in bdy_face_type)
    nperiodic +=1
    double1 = nor1*(xmax - xmin)
end
if ("periodicz" in bdy_face_type)
    nperiodic +=1
    double2 = nor2*(zmax - zmin)
end
if ("periodicy" in bdy_face_type)
    nperiodic +=1
    if (double1 == [0 ,0 ,0])
        double1 = nor3*(ymax-ymin)
    end
    if (double2 == [0 ,0 ,0])
        double2 = nor3*(ymax-ymin)
    end

end
if (nperiodic == 2) ## if we have triple periodicity this is handled above, for less than 2 periodic boundaries double periodicity does not apply
    ## find periodic planes
    plane1 = [0, 0, 0, 0]
    plane2 = [0, 0, 0, 0]
    plane1[1] = 1
    plane2_idx = 1
    plane1_idx = 2
    for i=2:8
        vec = [coords[1,1] - coords[i,1], coords[1,2] - coords[i,2], coords[1,3] - coords[i,3]]
        if (determine_colinearity(vec,double1) || determine_colinearity(vec,double2) || determine_colinearity(vec,abs.(double1+double2)))
            plane1[plane1_idx] = i
            plane1_idx += 1
        else
            plane2[plane2_idx] = i
            plane2_idx += 1
        end
    end
    target_idx = plane1[1]
    interval = [0, 0, 0]
    for i =1:3
        interval[i] = plane1[i+1]
    end
    #@info coords[target_idx], coords[target_idx,2],z[target_idx]
    for i in interval
        ip_kill=i
        ip_dest=target_idx
        for e=1:nelem
            for ii=1:ngl
                for jj=1:ngl
                    for kk=1:ngl
                        if (connijk[e,ii,jj,kk] == ip_kill)
                            connijk[e,ii,jj,kk] = ip_dest
                        end
                    end
                end
            end
        end
        for i=ip_kill:npoin-1
            coords[i,1] = coords[i+1,1]
            coords[i,2] = coords[i+1,2]
            coords[i,3] = coords[i+1,3]
        end
        npoin = npoin-1
        for iface =1:size(poin_in_bdy_face,1)
            for kk=1:ngl
                for ll=1:ngl
                    val = poin_in_bdy_face[iface,kk,ll]
                    if (val > ip_kill)
                        poin_in_bdy_face[iface,kk,ll] = val - 1
                    end
                    if (val == ip_kill)
                        poin_in_bdy_face[iface,kk,ll] = ip_dest
                    end
                end
            end
        end
        for e=1:nelem
            for ii=1:ngl
                for jj=1:ngl
                    for kk=1:ngl
                        ipp = connijk[e,ii,jj,kk]
                        if (ipp > ip_kill)
                            connijk[e,ii,jj,kk] -= 1
                        end
                    end
                end
            end
        end
        for i=1:size(interval,1)
            if (interval[i] > ip_kill)
                interval[i] -= 1
            elseif (interval[i] == ip_kill)
                interval[i] = 0
            end
        end
        for i=1:size(plane2,1)
            if (plane2[i] > ip_kill)
                plane2[i] -= 1
            elseif (plane2[i] == ip_kill)
                plane2[i] = ip_dest
            end
        end
    end
    target_idx = plane2[1]
    p2_idx = 1
    x_dest = coords[target_idx,1]
    y_dest = coords[target_idx,2]
    z_dest = coords[target_idx,3]
    for i=2:size(plane2,1)
        ii= plane2[i]
        x_dest = min(x_dest,coords[ii,1])
        y_dest = min(y_dest,coords[ii,2])
        z_dest = min(z_dest,coords[ii,3])
    end
    #make sure to pick a corner consistent with the vtk unwrap
    #=for i=2:size(plane2,1)
    ii = plane2[i]
    xt = coords[target_idx,1]
    yt = coords[target_idx,2]
    zt = coords[target_idx,3]
    xi = coords[ii,1]
    yi = coords[ii,2]
    zi = coords[ii,3]
    if (yi == 0 && yt == 0 && zi == 0 && zt == 0)
    comp1 = xi < xt
    elseif (yi == 0 && yt == 0)
    comp1 = xi*abs(zi) < xt*abs(zt)
    elseif (zi == 0 && zt == 0) 
    comp1 = xi*abs(zi) < xt*abs(yt)
    else
    comp1 = xi*abs(yi*zi) < xt*abs(yt*zt)
    end
    if (xt == xi && xi < 0.0) comp1 = !(comp1) end

    if (xi ==0 && xt == 0 && zi == 0 && zt ==0)
    comp2 = yi < yt
    elseif (xi == 0 && xt == 0)
    comp2 = yi*abs(zi) < yt*abs(zt)
    elseif (zi == 0 && zt == 0)
    comp2 = yi*abs(xi) < yt*abs(xt)
    else
    comp2 = yi*abs(xi*zi) < yt*abs(xt*zt)
    end
    if (yt == yi && yi < 0.0) comp2 = !(comp2) end 
    if (xi == 0 && xt == 0 && yi == 0 && yt ==0)
    comp3 = zi < zt
    elseif (xi == 0 && xt == 0)
    comp3 = zi*abs(yi) < zt*abs(yt)
    elseif (yi == 0 && yt == 0)
    comp3 = zi*abs(xi) < zt*abs(xt)
    else
    comp3 = zi*abs(xi*yi) < zt*abs(xt*yt)
    end
    if (zt == zi && zi < 0.0) comp3 = !(comp3) end
    xi,yi,zi,xt,yt,zt, comp1,comp2,comp3
    if (comp1 || comp2 || comp3)#(comp1 && comp2) || (comp1 && comp3) || (comp2 && comp3)
    @info xi,yi,zi,xt,yt,zt, comp1,comp2,comp3
    target_idx = plane2[i]
    p2_idx = i
    end
    end=#
    #@info x_dest, y_dest, z_dest

    interval = [0, 0, 0]
    for i =1:3
        interval[i] = plane2[i+1]
    end
    for i in interval
        ip_kill=i
        ip_dest=target_idx
        for e=1:nelem
            for ii=1:ngl
                for jj=1:ngl
                    for kk=1:ngl
                        if (connijk[e,ii,jj,kk] == ip_kill)
                            connijk[e,ii,jj,kk] = ip_dest
                        end
                    end
                end
            end
        end
        coords[ip_dest,1] = x_dest
        coords[ip_dest,2] = y_dest
        coords[ip_dest,3] = z_dest
        for i=ip_kill:npoin-1
            coords[i,1] = coords[i+1,1]
            coords[i,2] = coords[i+1,2]
            coords[i,3] = coords[i+1,3]
        end
        npoin = npoin-1
        for iface =1:size(poin_in_bdy_face,1)
            for kk=1:ngl
                for ll=1:ngl
                    val = poin_in_bdy_face[iface,kk,ll]
                    if (val > ip_kill)
                        poin_in_bdy_face[iface,kk,ll] = val - 1
                    end
                    if (val == ip_kill)
                        poin_in_bdy_face[iface,kk,ll] = ip_dest
                    end
                end
            end
        end
        for e=1:nelem
            for ii=1:ngl
                for jj=1:ngl
                    for kk=1:ngl
                        ipp = connijk[e,ii,jj,kk]
                        if (ipp > ip_kill)
                            connijk[e,ii,jj,kk] -= 1
                        end
                    end
                end
            end
        end
        for i=1:size(interval,1)
            if (interval[i] > ip_kill)
                interval[i] -= 1
            elseif (interval[i] == ip_kill)
                interval[i] = 0
            end
        end
        if (target_idx > ip_kill)
            target_idx -=1
        end
    end

end
### done with double periodicity
### do regular periodicity
per1_points = []
per2_points = []
per3_points = []
for iface_bdy =1:size(bdy_face_type,1)
    for k=1:ngl
        for l=1:ngl
            ip = poin_in_bdy_face[iface_bdy,k,l]
            if (bdy_face_type[iface_bdy] == "periodicx")
                per1_points = [per1_points; ip]
            elseif (bdy_face_type[iface_bdy] == "periodicz")
                per2_points = [per2_points; ip]
            elseif (bdy_face_type[iface_bdy] == "periodicy")
                per3_points = [per3_points; ip]
            end
        end
    end
end
### remove duplicates
unique!(per1_points)
unique!(per2_points)
unique!(per3_points)
### work single periodicity on the periodicx boundaries if applicable
if (size(per1_points,1) > 1) 
    for i=1:size(per1_points,1)
        found = false
        i1 = i+1
        while (i1 <= size(per1_points,1) && found == false)
            ip = per1_points[i]
            ip1 = per1_points[i1]
            vec = [coords[ip,1] - coords[ip1,1], coords[ip,2] - coords[ip1,2], coords[ip,3] - coords[ip1,3]]
            if (determine_colinearity(vec, nor1))
                #@info "found a match", vec, nor1, ip, ip1, coords[ip,1], coords[ip1,1]
                found = true
                xt = coords[ip1,1]
                yt = coords[ip1,2]
                zt = coords[ip1,3]
                xi = coords[ip,1]
                yi = coords[ip,2]
                zi = coords[ip,3]
                if (yi == 0 && yt == 0 && zi == 0 && zt == 0)
                    comp1 = xi < xt
                elseif (yi == 0 && yt == 0)
                    comp1 = xi*abs(zi) < xt*abs(zt)
                elseif (zi == 0 && zt == 0)
                    comp1 = xi*abs(zi) < xt*abs(yt)
                else
                    comp1 = xi*abs(yi*zi) < xt*abs(yt*zt)
                end
                if (xi ==0 && xt == 0 && zi == 0 && zt ==0)
                    comp2 = yi < yt
                elseif (xi == 0 && xt == 0)
                    comp2 = yi*abs(zi) < yt*abs(zt)
                elseif (zi == 0 && zt == 0)
                    comp2 = yi*abs(xi) < yt*abs(xt)
                else
                    comp2 = yi*abs(xi*zi) < yt*abs(xt*zt)
                end
                if (xi == 0 && xt == 0 && yi == 0 && yt ==0)
                    comp3 = zi < zt
                elseif (xi == 0 && xt == 0)
                    comp3 = zi*abs(yi) < zt*abs(yt)
                elseif (yi == 0 && yt == 0)
                    comp3 = zi*abs(xi) < zt*abs(xt)
                else
                    comp3 = zi*abs(xi*yi) < zt*abs(xt*yt)
                end
                #cond1 = coords[ip,1]*abs(coords[ip,2])*abs(coords[ip,3]) < coords[ip1,1]*abs(coords[ip1,2])*abs(coords[ip1,3])
                #cond2 = coords[ip,2]*abs(coords[ip,1])*abs(coords[ip,3]) < coords[ip1,2]*abs(coords[ip1,1])*abs(coords[ip1,3])
                #cond3 = coords[ip,3]*abs(coords[ip,1])*abs(coords[ip,2]) < coords[ip1,3]*abs(coords[ip1,1])*abs(coords[ip1,2])
                if (comp1 || comp2 || comp3)    
                    ip_dest = ip
                    ip_kill = ip1
                else
                    ip_dest = ip1
                    ip_kill = ip
                end
                connijk_spare .= connijk .- ip_kill
                poin_in_bdy_face_spare .= poin_in_bdy_face .- ip_kill
                for e=1:nelem
                    for ii=1:ngl
                        for jj=1:ngl
                            for kk=1:ngl
                                ipp = connijk_spare[e,ii,jj,kk]
                                if (ipp == 0)
                                    connijk[e,ii,jj,kk] = ip_dest
                                end
                            end
                        end
                    end
                end
                for iface = 1:size(poin_in_bdy_face,1)
                    for kk =1:ngl
                        for ll=1:ngl
                            ipp = poin_in_bdy_face_spare[iface,kk,ll]
                            if (ipp == 0)
                                poin_in_bdy_face[iface,kk,ll] = ip_dest
                            end
                        end
                    end
                end
                if !(ip_kill in connijk)
                    connijk_spare .= connijk .- ip_kill
                    poin_in_bdy_face_spare .= poin_in_bdy_face .- ip_kill
                    x_spare .= coords[:,1]
                    y_spare .= coords[:,2]
                    z_spare .= coords[:,3]
                    coords[ip_dest,1] = min(coords[ip_kill,1], coords[ip_dest,1])
                    coords[ip_dest,2] = min(coords[ip_kill,2], coords[ip_dest,2])
                    coords[ip_dest,3] = min(coords[ip_kill,3], coords[ip_dest,3])
                    @view(coords[ip_kill:npoin-1,1]) .= @view(x_spare[ip_kill+1:npoin])
                    @view(coords[ip_kill:npoin-1,2]) .= @view(y_spare[ip_kill+1:npoin])
                    @view(coords[ip_kill:npoin-1,3]) .= @view(z_spare[ip_kill+1:npoin])
                    x_spare .= coords[:,1]
                    y_spare .= coords[:,2]
                    z_spare .= coords[:,3]
                    npoin = npoin-1
                    mesh.npoin = mesh.npoin-1
                    for iface =1:size(poin_in_bdy_face,1)
                        for kk=1:ngl
                            for ll=1:ngl
                                val = poin_in_bdy_face_spare[iface,kk,ll]
                                if (val > 0)
                                    poin_in_bdy_face[iface,kk,ll] -= 1
                                end
                            end
                        end
                    end
                    for e=1:nelem
                        for ii=1:ngl
                            for jj=1:ngl
                                for kk=1:ngl
                                    ipp = connijk_spare[e,ii,jj,kk]
                                    comp = ipp > 0
                                    if (comp)
                                        connijk[e,ii,jj,kk] -= 1
                                    end
                                end
                            end
                        end
                    end
                    for ii=1:size(per1_points,1)
                        if (per1_points[ii] == ip_kill)
                            per1_points[ii] = ip_dest
                        end
                    end
                    for ii=1:size(per2_points,1)
                        if (per2_points[ii] == ip_kill)
                            per2_points[ii] = ip_dest
                        end
                    end
                    for ii=1:size(per3_points,1)
                        if (per3_points[ii] == ip_kill)
                            per3_points[ii] = ip_dest
                        end
                    end
                    for ii=1:size(per1_points,1)
                        if (per1_points[ii] > ip_kill)
                            per1_points[ii] -= 1
                        end
                    end
                    for ii=1:size(per2_points,1)
                        if (per2_points[ii] > ip_kill)
                            per2_points[ii] -= 1
                        end
                    end
                    for ii=1:size(per3_points,1)
                        if (per3_points[ii] > ip_kill)
                            per3_points[ii] -= 1
                        end
                    end
                end
            else
                i1 += 1
            end

        end
    end
end
unique!(per1_points)
unique!(per2_points)
unique!(per3_points)
### work on single periodicity on the periodicz boundaries if applicable
if (size(per2_points,1) > 1)    
    for i=1:size(per2_points,1)
        found = false
        i1 = i+1
        while (i1 <= size(per2_points,1) && found == false)
            ip = per2_points[i]
            ip1 = per2_points[i1]
            vec = [coords[ip,1] - coords[ip1,1], coords[ip,2] - coords[ip1,2], coords[ip,3] - coords[ip1,3]]
            if (determine_colinearity(vec, nor2))
                #@info "found a match", vec, nor2, ip, ip1,coords[ip,3], coords[ip1,3], coords[ip,2], coords[ip1,2], coords[ip,1], coords[ip1,1]
                found = true
                xt = coords[ip1,1]
                yt = coords[ip1,2]
                zt = coords[ip1,3]
                xi = coords[ip,1]
                yi = coords[ip,2]
                zi = coords[ip,3]
                if (yi == 0 && yt == 0 && zi == 0 && zt == 0)
                    comp1 = xi < xt
                elseif (yi == 0 && yt == 0)
                    comp1 = xi*abs(zi) < xt*abs(zt)
                elseif (zi == 0 && zt == 0)
                    comp1 = xi*abs(zi) < xt*abs(yt)
                else
                    comp1 = xi*abs(yi*zi) < xt*abs(yt*zt)
                end
                if (xi ==0 && xt == 0 && zi == 0 && zt ==0)
                    comp2 = yi < yt
                elseif (xi == 0 && xt == 0)
                    comp2 = yi*abs(zi) < yt*abs(zt)
                elseif (zi == 0 && zt == 0)
                    comp2 = yi*abs(xi) < yt*abs(xt)
                else
                    comp2 = yi*abs(xi*zi) < yt*abs(xt*zt)
                end
                if (xi == 0 && xt == 0 && yi == 0 && yt ==0)
                    comp3 = zi < zt
                elseif (xi == 0 && xt == 0)
                    comp3 = zi*abs(yi) < zt*abs(yt)
                elseif (yi == 0 && yt == 0)
                    comp3 = zi*abs(xi) < zt*abs(xt)
                else
                    comp3 = zi*abs(xi*yi) < zt*abs(xt*yt)
                end
                #cond1 = coords[ip,1]*abs(coords[ip,2])*abs(coords[ip,3]) < coords[ip1,1]*abs(coords[ip1,2])*abs(coords[ip1,3])
                #cond2 = coords[ip,2]*abs(coords[ip,1])*abs(coords[ip,3]) < coords[ip1,2]*abs(coords[ip1,1])*abs(coords[ip1,3])
                #cond3 = coords[ip,3]*abs(coords[ip,1])*abs(coords[ip,2]) < coords[ip1,3]*abs(coords[ip1,1])*abs(coords[ip1,2])
                if (comp1 || comp2 || comp3)
                    ip_dest = ip
                    ip_kill = ip1
                else
                    ip_dest = ip1
                    ip_kill = ip
                end
                connijk_spare .= connijk .- ip_kill
                poin_in_bdy_face_spare .= poin_in_bdy_face .- ip_kill
                for e=1:nelem
                    for ii=1:ngl
                        for jj=1:ngl
                            for kk=1:ngl
                                ipp = connijk_spare[e,ii,jj,kk]
                                if (ipp == 0)
                                    connijk[e,ii,jj,kk] = ip_dest
                                end
                            end
                        end
                    end
                end
                for iface = 1:size(poin_in_bdy_face,1)
                    for kk =1:ngl
                        for ll=1:ngl
                            ipp = poin_in_bdy_face_spare[iface,kk,ll]
                            if (ipp == 0)
                                poin_in_bdy_face[iface,kk,ll] = ip_dest
                            end
                        end
                    end
                end
                if !(ip_kill in connijk)
                    connijk_spare .= connijk .- ip_kill
                    poin_in_bdy_face_spare .= poin_in_bdy_face .- ip_kill
                    x_spare .= coords[:,1]
                    y_spare .= coords[:,2]
                    z_spare .= coords[:,3]
                    coords[ip_dest,1] = min(coords[ip_kill,1], coords[ip_dest,1])
                    coords[ip_dest,2] = min(coords[ip_kill,2], coords[ip_dest,2])
                    coords[ip_dest,3] = min(coords[ip_kill,3], coords[ip_dest,3])
                    @view(coords[ip_kill:npoin-1,1]) .= @view(x_spare[ip_kill+1:npoin])
                    @view(coords[ip_kill:npoin-1,2]) .= @view(y_spare[ip_kill+1:npoin])
                    @view(coords[ip_kill:npoin-1,3]) .= @view(z_spare[ip_kill+1:npoin])
                    x_spare .= coords[:,1]
                    y_spare .= coords[:,2]
                    z_spare .= coords[:,3]
                    npoin = npoin-1
                    for iface =1:size(poin_in_bdy_face,1)
                        for kk=1:ngl
                            for ll=1:ngl
                                val = poin_in_bdy_face_spare[iface,kk,ll]
                                if (val > 0)
                                    poin_in_bdy_face[iface,kk,ll] -= 1
                                end
                            end
                        end
                    end

                    for e=1:nelem
                        for ii=1:ngl
                            for jj=1:ngl
                                for kk=1:ngl
                                    ipp = connijk_spare[e,ii,jj,kk]
                                    if (ipp > 0)
                                        connijk[e,ii,jj,kk] -= 1
                                    end
                                end
                            end
                        end
                    end
                    for ii=1:size(per1_points,1)
                        if (per1_points[ii] == ip_kill)
                            per1_points[ii] = ip_dest
                        end
                    end
                    for ii=1:size(per2_points,1)
                        if (per2_points[ii] == ip_kill)
                            per2_points[ii] = ip_dest
                        end
                    end
                    for ii=1:size(per3_points,1)
                        if (per3_points[ii] == ip_kill)
                            per3_points[ii] = ip_dest
                        end
                    end
                    for ii=1:size(per1_points,1)
                        if (per1_points[ii] > ip_kill)
                            per1_points[ii] -= 1
                        end
                    end
                    for ii=1:size(per2_points,1)
                        if (per2_points[ii] > ip_kill)
                            per2_points[ii] -= 1
                        end
                    end
                    for ii=1:size(per3_points,1)
                        if (per3_points[ii] > ip_kill)
                            per3_points[ii] -= 1
                        end
                    end
                end
            else
                i1 += 1
            end

        end
    end
end
unique!(per1_points)
unique!(per2_points)
unique!(per3_points)
### work on single periodicity on the periodicy boundaries if applicable
if (size(per3_points,1) > 1)
    for i=1:size(per3_points,1)
        found = false
        i1 = i+1
        while (i1 <= size(per3_points,1) && found == false)
            ip = per3_points[i]
            ip1 = per3_points[i1]
            vec = [coords[ip,1] - coords[ip1,1], coords[ip,2] - coords[ip1,2], coords[ip,3] - coords[ip1,3]]
            if (determine_colinearity(vec, nor3))
                #@info "found a match", vec, nor3, ip, ip1, coords[ip,2], coords[ip1,2], coords[ip,1], coords[ip1,1], coords[ip,3], coords[ip1,3]
                found = true
                xt = coords[ip1,1]
                yt = coords[ip1,2]
                zt = coords[ip1,3]
                xi = coords[ip,1]
                yi = coords[ip,2]
                zi = coords[ip,3]
                if (yi == 0 && yt == 0 && zi == 0 && zt == 0)
                    comp1 = xi < xt
                elseif (yi == 0 && yt == 0)
                    comp1 = xi*abs(zi) < xt*abs(zt)
                elseif (zi == 0 && zt == 0)
                    comp1 = xi*abs(zi) < xt*abs(yt)
                else
                    comp1 = xi*abs(yi*zi) < xt*abs(yt*zt)
                end
                if (xi ==0 && xt == 0 && zi == 0 && zt ==0)
                    comp2 = yi < yt
                elseif (xi == 0 && xt == 0)
                    comp2 = yi*abs(zi) < yt*abs(zt)
                elseif (zi == 0 && zt == 0)
                    comp2 = yi*abs(xi) < yt*abs(xt)
                else
                    comp2 = yi*abs(xi*zi) < yt*abs(xt*zt)
                end
                if (xi == 0 && xt == 0 && yi == 0 && yt ==0)
                    comp3 = zi < zt
                elseif (xi == 0 && xt == 0)
                    comp3 = zi*abs(yi) < zt*abs(yt)
                elseif (yi == 0 && yt == 0)
                    comp3 = zi*abs(xi) < zt*abs(xt)
                else
                    comp3 = zi*abs(xi*yi) < zt*abs(xt*yt)
                end
                #cond1 = coords[ip,1]*abs(coords[ip,2])*abs(coords[ip,3]) < coords[ip1,1]*abs(coords[ip1,2])*abs(coords[ip1,3])
                #cond2 = coords[ip,2]*abs(coords[ip,1])*abs(coords[ip,3]) < coords[ip1,2]*abs(coords[ip1,1])*abs(coords[ip1,3])
                #cond3 = coords[ip,3]*abs(coords[ip,1])*abs(coords[ip,2]) < coords[ip1,3]*abs(coords[ip1,1])*abs(coords[ip1,2])
                if (comp1 || comp2 || comp3)
                    ip_dest = ip
                    ip_kill = ip1
                else
                    ip_dest = ip1
                    ip_kill = ip
                end
                connijk_spare .= connijk .- ip_kill
                poin_in_bdy_face_spare .= poin_in_bdy_face .- ip_kill
                for e=1:nelem
                    for ii=1:ngl
                        for jj=1:ngl
                            for kk=1:ngl
                                ipp = connijk_spare[e,ii,jj,kk]
                                if (ipp == 0)
                                    connijk[e,ii,jj,kk] = ip_dest
                                end
                            end
                        end
                    end
                end
                for iface = 1:size(poin_in_bdy_face,1)
                    for kk =1:ngl
                        for ll=1:ngl
                            ipp = poin_in_bdy_face_spare[iface,kk,ll]
                            if (ipp == 0)
                                poin_in_bdy_face[iface,kk,ll] = ip_dest
                            end
                        end
                    end
                end
                if !(ip_kill in connijk)
                    connijk_spare .= connijk .- ip_kill
                    poin_in_bdy_face_spare .= poin_in_bdy_face .- ip_kill
                    x_spare .= coords[:,1]
                    y_spare .= coords[:,2]
                    z_spare .= coords[:,3]
                    coords[ip_dest,1] = min(coords[ip_kill,1], coords[ip_dest,1])
                    coords[ip_dest,2] = min(coords[ip_kill,2], coords[ip_dest,2])
                    coords[ip_dest,3] = min(coords[ip_kill,3], coords[ip_dest,3])
                    @view(coords[ip_kill:npoin-1,1]) .= @view(x_spare[ip_kill+1:npoin])
                    @view(coords[ip_kill:npoin-1,2]) .= @view(y_spare[ip_kill+1:npoin])
                    @view(coords[ip_kill:npoin-1,3]) .= @view(z_spare[ip_kill+1:npoin])
                    x_spare .= coords[:,1]
                    y_spare .= coords[:,2]
                    z_spare .= coords[:,3]
                    npoin = npoin-1
                    #=coords[ip_dest] = min(coords[ip_kill,1],coords[ip_dest,1])
                    coords[ip_dest,2] = min(coords[ip_kill,2],coords[ip_dest,2])
                    coords[ip_dest,3] = min(coords[ip_kill,3],coords[ip_dest,3])
                    for ii=ip_kill:npoin-1
                    coords[ii,1] = coords[ii+1,1]
                    coords[ii,2] = coords[ii+1,2]
                    coords[ii,3] = coords[ii+1,3]
                    end
                    npoin = npoin-1=#
                    for iface =1:size(poin_in_bdy_face,1)
                        for kk=1:ngl
                            for ll=1:ngl
                                val = poin_in_bdy_face_spare[iface,kk,ll]
                                if (val > 0)
                                    poin_in_bdy_face[iface,kk,ll] -= 1
                                end
                            end
                        end
                    end

                    for e=1:nelem
                        for ii=1:ngl
                            for jj=1:ngl
                                for kk=1:ngl
                                    ipp = connijk_spare[e,ii,jj,kk]
                                    if (ipp > 0)
                                        connijk[e,ii,jj,kk] -= 1
                                    end
                                end
                            end
                        end
                    end
                    for ii=1:size(per1_points,1)
                        if (per1_points[ii] == ip_kill)
                            per1_points[ii] = ip_dest
                        end
                    end
                    for ii=1:size(per2_points,1)
                        if (per2_points[ii] == ip_kill)
                            per2_points[ii] = ip_dest
                        end
                    end
                    for ii=1:size(per3_points,1)
                        if (per3_points[ii] == ip_kill)
                            per3_points[ii] = ip_dest
                        end
                    end
                    for ii=1:size(per1_points,1)
                        if (per1_points[ii] > ip_kill)
                            per1_points[ii] -= 1
                        end
                    end
                    for ii=1:size(per2_points,1)
                        if (per2_points[ii] > ip_kill)
                            per2_points[ii] -= 1
                        end
                    end
                    for ii=1:size(per3_points,1)
                        if (per3_points[ii] > ip_kill)
                            per3_points[ii] -= 1
                        end
                    end
                end
            else
                i1 += 1
            end

        end
    end
end
unique!(per1_points)
unique!(per2_points)
unique!(per3_points)
for e=1:nelem
    for i=1:ngl
        for j=1:ngl
            for k=1:ngl
                ip = connijk[e,i,j,k]
            end
        end
    end
end
mesh.npoin = npoin
else
#
# 1D periodicity
#
if inputs[:AD] != FD()
    ip_dest = 1
    ip_kill = npoin_linear
    for e=1:nelem
        for i=1:ngl
            if (connijk[e,i,1] == ip_kill)
                connijk[e,i,1] = ip_dest
            elseif (connijk[e,i,1] > ip_kill)
                connijk[e,i,1] -= 1
            end
        end
    end
    for ip=ip_kill:npoin-1
        coords[ip,1] = coords[ip+1,1]
    end
    npoin = npoin-1
end
end
mesh.npoin = npoin
if mesh.rank == 0
    @info " periodicity_restructure!"
end
end

function restructure4periodicity_2D(mesh, norm, periodic_direction)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)
    
    per_ip = Int[]
    ngl = mesh.ngl
    for iedge_bdy =1:size(mesh.bdy_edge_type,1)
        for k=1:ngl
            ip = mesh.poin_in_bdy_edge[iedge_bdy,k]
            if (mesh.bdy_edge_type[iedge_bdy] == periodic_direction)
                per_ip = [per_ip; ip]
            end
        end
    end

    if (mesh.lLaguerre)
        e_iter = 1
        for iedge_bdy =1:size(mesh.bdy_edge_type,1)
            if (mesh.bdy_edge_type[iedge_bdy] == "Laguerre")
                if (mesh.poin_in_bdy_edge[iedge_bdy,1] in per_ip)
                    for k = 2:mesh.ngr
                        ip = mesh.connijk_lag[e_iter,1,k]
                        per_ip = [per_ip; ip]
                    end
                elseif (mesh.poin_in_bdy_edge[iedge_bdy,mesh.ngl] in per_ip)
                    for k = 2:mesh.ngr
                        ip = mesh.connijk_lag[e_iter,mesh.ngl,k]
                        per_ip = [per_ip; ip]
                    end
                end
                e_iter += 1
            end
        end
    end

    ### remove duplicates
    unique!(per_ip)
    x_local  = mesh.x[per_ip]
    y_local  = mesh.y[per_ip]
    per_gip  = mesh.ip2gip[per_ip]
    ip_owner = mesh.gip2owner[per_ip]
    # @info  mesh.x[per_ip]

    # Gather arrays onto the root processor (rank 0)
    root = 0
	# Gather per_gip
    buffer_sz::Int32    = size(per_ip, 1)
    # @info rank, buffer_sz
    recv_counts  = MPI.Gather(buffer_sz, 0, comm)
    # total_counts = sum(recv_counts)
    # if total_counts == 0
        # return
    # end
    # else
    # if total_counts == 0
    #     return
    # end

    x_gather     = MPI.gather(x_local, comm)
    y_gather     = MPI.gather(y_local, comm)
    gathered_per = MPI.gather(per_gip, comm)
    owner_gather = MPI.gather(ip_owner, comm)
    if mesh.rank == root

    # On the root processor, combine and remove duplicates
        # Concatenate gathered arrays
        x              = vcat(x_gather...)
        y              = vcat(y_gather...)
        global_per_gip = vcat(gathered_per...)
        owner          = vcat(owner_gather...)

        sz = size(global_per_gip,1)
		for i = 1:sz
            i1 = i+1
            for i1 = (i+1):sz
                if global_per_gip[i] == global_per_gip[i1]
                    continue
                end
                vec = [x[i] - x[i1], y[i] - y[i1]]
                # @info vec, norm
                if (determine_colinearity(vec, norm))
                    xt = x[i1]
                    yt = y[i1]
                    xi = x[i]
                    yi = y[i]
                    if (yi == 0 && yt == 0)
                        comp1 = xi < xt
                    else
                        comp1 = xi*abs(yi) < xt*abs(yt)
                    end
                    if (xi ==0 && xt == 0)
                        comp2 = yi < yt
                    else
                        comp2 = yi*abs(xi) < yt*abs(xt)
                    end
                    # @info "found", global_per_gip[i], global_per_gip[i1]
                    if (comp1 || comp2)
                        global_per_gip[i1] = global_per_gip[i]
                        if owner[i1] != owner[i]
                            owner[i1] = owner[i]
                        end
                    else
                        global_per_gip[i] = global_per_gip[i1]
                        if owner[i1] != owner[i]
                            owner[i] = owner[i1]
                        end
                    end
                    # break
                else
                    continue
                end
            end
        end
		# do something for global_per_gip
        s_gip_vbuf   = VBuffer(global_per_gip, recv_counts)
        s_owner_vbuf = VBuffer(owner, recv_counts)
    else
        s_gip_vbuf   = VBuffer(nothing)
        s_owner_vbuf = VBuffer(nothing)
    end
    MPI.Barrier(comm)
    per_ip_updated = MPI.Scatterv!(s_gip_vbuf,zeros(eltype(per_gip), buffer_sz), 0, comm)
    owner_updated  = MPI.Scatterv!(s_owner_vbuf,zeros(eltype(ip_owner), buffer_sz), 0, comm)
    # per_ip_updated = MPI.Scatterv!(global_per_gip,buffer_sz, 0, comm)

    mesh.ip2gip[per_ip]    .= per_ip_updated
    mesh.gip2owner[per_ip] .= owner_updated
end


function restructure4periodicity_3D_sorted!(mesh, norm, periodic_direction)
    comm    = MPI.COMM_WORLD
    rank    = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)
    ngl     = mesh.ngl
    ngl2    = ngl * ngl

    # Collect periodic boundary nodes, tagging each as a face-corner node
    # (both k and l at an extreme index) or not.
    # Child element corners at AMR refinement midpoints are corner-type;
    # the coarse parent's mid-edge GLL node at the same position is non-corner-type.
    # Keeping them in separate groups prevents the wrong periodic merge.
    ip_is_corner = Dict{Int, Bool}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        mesh.bdy_face_type[iface_bdy] == periodic_direction || continue
        for k = 1:ngl, l = 1:ngl
            ip = mesh.poin_in_bdy_face[iface_bdy, k, l]
            corner = (k == 1 || k == ngl) && (l == 1 || l == ngl)
            ip_is_corner[ip] = get(ip_is_corner, ip, false) || corner
        end
    end
    per_ip = collect(keys(ip_is_corner))



    x_local  = mesh.x[per_ip]
    y_local  = mesh.y[per_ip]
    z_local  = mesh.z[per_ip]
    per_gip  = mesh.ip2gip[per_ip]
    ip_owner = mesh.gip2owner[per_ip]
    is_corner_local = Int8[ip_is_corner[ip] ? 1 : 0 for ip in per_ip]

    # Gather arrays onto the root processor (rank 0)
    root = 0

    buffer_sz::Int32    = size(per_ip, 1)
    recv_counts  = MPI.Gather(buffer_sz, 0, comm)

    x_gather      = MPI.gather(x_local, comm)
    y_gather      = MPI.gather(y_local, comm)
    z_gather      = MPI.gather(z_local, comm)
    gathered_per  = MPI.gather(per_gip, comm)
    owner_gather  = MPI.gather(ip_owner, comm)
    corner_gather = MPI.gather(is_corner_local, comm)



    if mesh.rank == root

        # On the root processor, combine and remove duplicates
        # Concatenate gathered arrays
        x              = vcat(x_gather...)
        y              = vcat(y_gather...)
        z              = vcat(z_gather...)
        global_per_gip  = vcat(gathered_per...)
        owner           = vcat(owner_gather...)
        global_is_corner = vcat(corner_gather...)
        # Build per-gip corner flag: a gip is "corner" if any contributing rank
        # saw it at a face-corner position (k and l both extreme).
        gip_is_corner = Dict{Int, Bool}()
        for (i, g) in enumerate(global_per_gip)
            gip_is_corner[g] = get(gip_is_corner, g, false) || (global_is_corner[i] == 1)
        end
        # ─────────────────────────────────────────────────────────────────────
        # Deduplicate by (x, y, z, is_corner) so that at NCF periodic positions a
        # child-corner node and a coarse-neighbor mid-edge node at the same coordinate
        # are both kept and paired independently with their correct periodic partners.
        coords = collect(zip(round.(x; digits=5), round.(y; digits=5), round.(z; digits=5), global_is_corner))
        uniq_idx = unique(i -> coords[i], eachindex(coords))
        un_gathered_x = collect(@view x[uniq_idx])
        un_gathered_y = collect(@view y[uniq_idx])
        un_gathered_z = collect(@view z[uniq_idx])
        un_updated_global_per_gip = collect(@view global_per_gip[uniq_idx])
        un_updated_global_owner   = collect(@view owner[uniq_idx])
        un_is_corner              = collect(@view global_is_corner[uniq_idx])


        changes_ip    = Dict{Int, Int}()
        changes_owner = Dict{Int, Int}()
        sz = length(uniq_idx)
        if sz > 0
            # Choose periodic axis and transverse axes.
            if periodic_direction == "periodicx"
                ax_arr = un_gathered_x; t1_arr = un_gathered_y; t2_arr = un_gathered_z
            elseif periodic_direction == "periodicy"
                ax_arr = un_gathered_y; t1_arr = un_gathered_x; t2_arr = un_gathered_z
            else
                ax_arr = un_gathered_z; t1_arr = un_gathered_x; t2_arr = un_gathered_y
            end
            axmin, axmax = extrema(ax_arr)

            # Build lookup for max-side nodes: (t1_rounded, t2_rounded, is_corner) → (gip, owner, idx).
            # Including is_corner in the key means child-corner nodes only pair with child-corner
            # nodes and coarse mid-edge nodes only pair with coarse mid-edge nodes.
            # Nodes with no matching counterpart of the same corner type (hanging nodes, NCF-periodic)
            # are left unmerged via the haskey miss below.
            max_node_lookup = Dict{NTuple{3,Float64}, NTuple{3,Int}}()
            for idx in eachindex(ax_arr)
                abs(ax_arr[idx] - axmax) < 1e-4 || continue
                key = (round(t1_arr[idx]; digits=5), round(t2_arr[idx]; digits=5),
                       Float64(un_is_corner[idx]))
                max_node_lookup[key] = (un_updated_global_per_gip[idx],
                                        un_updated_global_owner[idx], idx)
            end

            vec = fill!(similar(norm), 0.0)
            for idx_i in eachindex(ax_arr)
                abs(ax_arr[idx_i] - axmin) < 1e-4 || continue
                key = (round(t1_arr[idx_i]; digits=5), round(t2_arr[idx_i]; digits=5),
                       Float64(un_is_corner[idx_i]))
                # No max-side node at same position with same corner type →
                # NCF-periodic or hanging node; leave unmerged.
                haskey(max_node_lookup, key) || continue
                gip_j_val, _, idx_j = max_node_lookup[key]
                vec[1] = un_gathered_x[idx_i] - un_gathered_x[idx_j]
                vec[2] = un_gathered_y[idx_i] - un_gathered_y[idx_j]
                vec[3] = un_gathered_z[idx_i] - un_gathered_z[idx_j]
                gip_i = un_updated_global_per_gip[idx_i]
                gip_j = gip_j_val
                if (get(changes_ip, gip_i, gip_i) == gip_j) || get(changes_ip, gip_j, gip_j) == gip_i
                    continue
                end
                if (determine_colinearity(vec, norm))
                    xt = un_gathered_x[idx_j]
                    yt = un_gathered_y[idx_j]
                    zt = un_gathered_z[idx_j]
                    xi = un_gathered_x[idx_i]
                    yi = un_gathered_y[idx_i]
                    zi = un_gathered_z[idx_i]
                    if (yi == 0 && yt == 0 && zi == 0 && zt == 0)
                        comp1 = xi < xt
                    elseif (yi == 0 && yt == 0)
                        comp1 = xi*abs(zi) < xt*abs(zt)
                    elseif (zi == 0 && zt == 0)
                        comp1 = xi*abs(zi) < xt*abs(yt)
                    else
                        comp1 = xi*abs(yi*zi) < xt*abs(yt*zt)
                    end
                    if (xi ==0 && xt == 0 && zi == 0 && zt ==0)
                        comp2 = yi < yt
                    elseif (xi == 0 && xt == 0)
                        comp2 = yi*abs(zi) < yt*abs(zt)
                    elseif (zi == 0 && zt == 0)
                        comp2 = yi*abs(xi) < yt*abs(xt)
                    else
                        comp2 = yi*abs(xi*zi) < yt*abs(xt*zt)
                    end
                    if (xi == 0 && xt == 0 && yi == 0 && yt ==0)
                        comp3 = zi < zt
                    elseif (xi == 0 && xt == 0)
                        comp3 = zi*abs(yi) < zt*abs(yt)
                    elseif (yi == 0 && yt == 0)
                        comp3 = zi*abs(xi) < zt*abs(xt)
                    else
                        comp3 = zi*abs(xi*yi) < zt*abs(xt*yt)
                    end
                    if comp1 || comp2 || comp3
                        # j is the slave, i is the master
                        changes_ip[gip_j] = gip_i
                        if un_updated_global_owner[idx_j] != un_updated_global_owner[idx_i]
                            changes_owner[gip_j] = un_updated_global_owner[idx_i]
                            changes_owner[gip_i] = un_updated_global_owner[idx_i]
                        end
                    else
                        # i is the slave, j is the master
                        changes_ip[gip_i] = gip_j
                        if un_updated_global_owner[idx_j] != un_updated_global_owner[idx_i]
                            changes_owner[gip_i] = un_updated_global_owner[idx_j]
                            changes_owner[gip_j] = un_updated_global_owner[idx_j]
                        end
                    end
                end
                # If determine_colinearity fails for a position-matched pair it means
                # numerical noise; skip rather than crash.
            end
        end
        updated_global_per_gip = [get(changes_ip, x, x) for x in global_per_gip]
        updated_owner = [get(changes_owner, x, owner[i])  for (i, x) in enumerate(global_per_gip)]
        s_gip_vbuf   = VBuffer(updated_global_per_gip, recv_counts)
        s_owner_vbuf = VBuffer(updated_owner, recv_counts)
    else
        s_gip_vbuf   = VBuffer(nothing)
        s_owner_vbuf = VBuffer(nothing)
    end
    MPI.Barrier(comm)
    per_ip_updated = MPI.Scatterv!(s_gip_vbuf,zeros(eltype(per_gip), buffer_sz), 0, comm)
    owner_updated  = MPI.Scatterv!(s_owner_vbuf,zeros(eltype(ip_owner), buffer_sz), 0, comm)

    mesh.ip2gip[per_ip]    .= per_ip_updated
    mesh.gip2owner[per_ip] .= owner_updated
end


function restructure_el2gel_for_periodicity_3D!(mesh, _norm, periodic_direction)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    ngl  = mesh.ngl

    # Collect local periodic boundary face data: one entry per face
    per_faces_iel  = Int[]
    per_faces_gel  = eltype(mesh.el2gel)[]
    per_faces_own  = eltype(mesh.gel2owner)[]
    cx_local       = Float64[]
    cy_local       = Float64[]
    cz_local       = Float64[]

    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        if mesh.bdy_face_type[iface_bdy] == periodic_direction
            iel = mesh.bdy_face_in_elem[iface_bdy]
            gel = mesh.el2gel[iel]
            own = mesh.gel2owner[iel]
            cx = 0.0; cy = 0.0; cz = 0.0
            for k = 1:ngl, l = 1:ngl
                ip  = mesh.poin_in_bdy_face[iface_bdy, k, l]
                cx += mesh.x[ip]
                cy += mesh.y[ip]
                cz += mesh.z[ip]
            end
            npts = Float64(ngl * ngl)
            push!(per_faces_iel, iel)
            push!(per_faces_gel, gel)
            push!(per_faces_own, own)
            push!(cx_local, cx / npts)
            push!(cy_local, cy / npts)
            push!(cz_local, cz / npts)
        end
    end

    buffer_sz::Int32 = length(per_faces_iel)
    recv_counts  = MPI.Gather(buffer_sz, 0, comm)
    gel_gather   = MPI.gather(per_faces_gel, comm)
    own_gather   = MPI.gather(per_faces_own, comm)
    cx_gather    = MPI.gather(cx_local, comm)
    cy_gather    = MPI.gather(cy_local, comm)
    cz_gather    = MPI.gather(cz_local, comm)

    if rank == 0
        all_gel = vcat(gel_gather...)
        all_own = vcat(own_gather...)
        all_cx  = vcat(cx_gather...)
        all_cy  = vcat(cy_gather...)
        all_cz  = vcat(cz_gather...)

        changes_gel = Dict{eltype(mesh.el2gel),   eltype(mesh.el2gel)}()
        changes_own = Dict{eltype(mesh.el2gel),   eltype(mesh.gel2owner)}()

        if length(all_gel) > 0
            # Choose the periodic coordinate and transverse coordinates
            if periodic_direction == "periodicx"
                x1, x2, x3 = all_cy, all_cz, all_cx
            elseif periodic_direction == "periodicy"
                x1, x2, x3 = all_cx, all_cz, all_cy
            elseif periodic_direction == "periodicz"
                x1, x2, x3 = all_cx, all_cy, all_cz
            end

            x3min, x3max = extrema(x3)
            idx_min = findall(v -> AlmostEqual(v, x3min), x3)
            idx_max = findall(v -> AlmostEqual(v, x3max), x3)

            # Build centroid lookup for max-side faces: rounded (x1,x2) → (gel, own, idx)
            # Only conforming pairs (same centroid) are merged; NCF pairs are skipped.
            max_face_centroid = Dict{NTuple{2,Float64},Tuple{Int,Int,Int}}()
            for k in idx_max
                key = (round(x1[k]; digits=5), round(x2[k]; digits=5))
                max_face_centroid[key] = (all_gel[k], all_own[k], k)
            end

            for gi in idx_min
                key = (round(x1[gi]; digits=5), round(x2[gi]; digits=5))
                haskey(max_face_centroid, key) || continue  # NCF pair; no el2gel merge
                gel_j_v, own_j_v, gj = max_face_centroid[key]
                gel_i = get(changes_gel, all_gel[gi], all_gel[gi])
                gel_j = get(changes_gel, all_gel[gj], all_gel[gj])
                own_i = get(changes_own, all_gel[gi], all_own[gi])
                own_j = get(changes_own, all_gel[gj], all_own[gj])
                if gel_i <= gel_j
                    changes_gel[all_gel[gi]] = gel_i
                    changes_gel[all_gel[gj]] = gel_i
                    changes_own[all_gel[gi]] = own_i
                    changes_own[all_gel[gj]] = own_i
                else
                    changes_gel[all_gel[gi]] = gel_j
                    changes_gel[all_gel[gj]] = gel_j
                    changes_own[all_gel[gi]] = own_j
                    changes_own[all_gel[gj]] = own_j
                end
            end
        end

        updated_gel = [get(changes_gel, g, g) for g in all_gel]
        updated_own = [get(changes_own, g, all_own[i]) for (i, g) in enumerate(all_gel)]
        s_gel_vbuf  = VBuffer(updated_gel, recv_counts)
        s_own_vbuf  = VBuffer(updated_own, recv_counts)
    else
        s_gel_vbuf = VBuffer(nothing)
        s_own_vbuf = VBuffer(nothing)
    end

    MPI.Barrier(comm)
    gel_updated = MPI.Scatterv!(s_gel_vbuf, zeros(eltype(mesh.el2gel),   buffer_sz), 0, comm)
    own_updated = MPI.Scatterv!(s_own_vbuf, zeros(eltype(mesh.gel2owner), buffer_sz), 0, comm)

    for (i, iel) in enumerate(per_faces_iel)
        mesh.el2gel[iel]   = gel_updated[i]
        mesh.gel2owner[iel] = own_updated[i]
    end
end


function restructure_el2gel_for_periodicity_2D!(mesh, _norm, periodic_direction)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    ngl  = mesh.ngl

    # Collect local periodic boundary edge data: one entry per edge
    per_edges_iel = Int[]
    per_edges_gel = eltype(mesh.el2gel)[]
    cx_local      = Float64[]
    cy_local      = Float64[]

    for iedge_bdy = 1:size(mesh.bdy_edge_type, 1)
        if mesh.bdy_edge_type[iedge_bdy] == periodic_direction
            iel = mesh.bdy_edge_in_elem[iedge_bdy]
            gel = mesh.el2gel[iel]
            cx = 0.0; cy = 0.0
            for k = 1:ngl
                ip  = mesh.poin_in_bdy_edge[iedge_bdy, k]
                cx += mesh.x[ip]
                cy += mesh.y[ip]
            end
            push!(per_edges_iel, iel)
            push!(per_edges_gel, gel)
            push!(cx_local, cx / Float64(ngl))
            push!(cy_local, cy / Float64(ngl))
        end
    end

    buffer_sz::Int32 = length(per_edges_iel)
    recv_counts = MPI.Gather(buffer_sz, 0, comm)
    gel_gather  = MPI.gather(per_edges_gel, comm)
    cx_gather   = MPI.gather(cx_local, comm)
    cy_gather   = MPI.gather(cy_local, comm)

    if rank == 0
        all_gel = vcat(gel_gather...)
        all_cx  = vcat(cx_gather...)
        all_cy  = vcat(cy_gather...)

        changes_gel = Dict{eltype(mesh.el2gel), eltype(mesh.el2gel)}()

        if length(all_gel) > 0
            # Choose the periodic coordinate and transverse coordinate
            if periodic_direction == "periodicx"
                x1, x3 = all_cy, all_cx
            else  # periodicy / periodicz
                x1, x3 = all_cx, all_cy
            end

            x3min, x3max = extrema(x3)
            idx_min = findall(v -> AlmostEqual(v, x3min), x3)
            idx_max = findall(v -> AlmostEqual(v, x3max), x3)

            sort_min = sortperm(x1[idx_min])
            sort_max = sortperm(x1[idx_max])
            idx_min_sorted = idx_min[sort_min]
            idx_max_sorted = idx_max[sort_max]

            for i = 1:min(length(idx_min_sorted), length(idx_max_sorted))
                gel_i = all_gel[idx_min_sorted[i]]
                gel_j = all_gel[idx_max_sorted[i]]
                gel_i = get(changes_gel, gel_i, gel_i)
                gel_j = get(changes_gel, gel_j, gel_j)
                gel_min = min(gel_i, gel_j)
                changes_gel[gel_i] = gel_min
                changes_gel[gel_j] = gel_min
            end
        end

        updated_gel = [get(changes_gel, g, g) for g in all_gel]
        s_gel_vbuf  = VBuffer(updated_gel, recv_counts)
    else
        s_gel_vbuf = VBuffer(nothing)
    end

    MPI.Barrier(comm)
    gel_updated = MPI.Scatterv!(s_gel_vbuf, zeros(eltype(mesh.el2gel), buffer_sz), 0, comm)

    for (i, iel) in enumerate(per_edges_iel)
        mesh.el2gel[iel] = gel_updated[i]
    end
end


# ---------------------------------------------------------------------------
# Detect nonconforming periodic face pairs in 3D and MERGE them directly into
# the existing non_conforming_facets / ghost lists so that the standard
# DSS_nc_gather_rhs! / DSS_nc_scatter_rhs! handle them without any new code.
#
# Entry format is [ciel, piel, cfid, h1, h2] — identical to interior NCF:
#   cfid = child's face ID (1-6); parent's face is always opposite(cfid).
#
# Call collect_periodic_ncf_pairs_3D! once per direction, then call
# extend_ncf_ip_lists_3D! once to append new columns to IPc_list / IPp_list.
#
# Face-ID convention:
#   1 front  (j=ngl)   2 back   (j=1)
#   3 bottom (k=ngl)   4 top    (k=1)
#   5 right  (i=1)     6 left   (i=ngl)
#
# For periodicx: child on min-side → cfid=5, child on max-side → cfid=6
# For periodicy: child on min-side → cfid=2, child on max-side → cfid=1
# For periodicz: child on min-side → cfid=4, child on max-side → cfid=3
# ---------------------------------------------------------------------------
function collect_periodic_ncf_pairs_3D!(mesh, periodic_direction, elm2pelm)
    comm  = MPI.COMM_WORLD
    rank  = MPI.Comm_rank(comm)
    ngl   = mesh.ngl
    ngl2  = ngl * ngl

    # Determine face IDs for each side of this direction.
    if periodic_direction == "periodicx"
        fid_min = 5; fid_max = 6
    elseif periodic_direction == "periodicy"
        fid_min = 2; fid_max = 1
    else  # periodicz
        fid_min = 4; fid_max = 3
    end

    # -----------------------------------------------------------------
    # Phase 1 – each rank collects its periodic boundary faces.
    # -----------------------------------------------------------------
    local_iel  = Int[]
    local_gel  = eltype(elm2pelm)[]
    local_rk   = Int[]
    cx1_loc    = Float64[]
    cx2_loc    = Float64[]
    cx3_loc    = Float64[]
    x1mn_loc   = Float64[]
    x1mx_loc   = Float64[]
    x2mn_loc   = Float64[]
    x2mx_loc   = Float64[]

    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        mesh.bdy_face_type[iface_bdy] == periodic_direction || continue
        iel = mesh.bdy_face_in_elem[iface_bdy]
        gel = elm2pelm[iel]
        cx = 0.0; cy = 0.0; cz = 0.0
        x1mn = Inf; x1mx = -Inf; x2mn = Inf; x2mx = -Inf
        for kk = 1:ngl, ll = 1:ngl
            ip = mesh.poin_in_bdy_face[iface_bdy, kk, ll]
            xx = mesh.x[ip]; yy = mesh.y[ip]; zz = mesh.z[ip]
            cx += xx; cy += yy; cz += zz
            if periodic_direction == "periodicx"
                x1mn = min(x1mn, yy); x1mx = max(x1mx, yy)
                x2mn = min(x2mn, zz); x2mx = max(x2mx, zz)
            elseif periodic_direction == "periodicy"
                x1mn = min(x1mn, xx); x1mx = max(x1mx, xx)
                x2mn = min(x2mn, zz); x2mx = max(x2mx, zz)
            else
                x1mn = min(x1mn, xx); x1mx = max(x1mx, xx)
                x2mn = min(x2mn, yy); x2mx = max(x2mx, yy)
            end
        end
        npts = Float64(ngl2)
        push!(local_iel, iel); push!(local_gel, gel); push!(local_rk, rank)
        if periodic_direction == "periodicx"
            push!(cx1_loc, cy/npts); push!(cx2_loc, cz/npts); push!(cx3_loc, cx/npts)
        elseif periodic_direction == "periodicy"
            push!(cx1_loc, cx/npts); push!(cx2_loc, cz/npts); push!(cx3_loc, cy/npts)
        else
            push!(cx1_loc, cx/npts); push!(cx2_loc, cy/npts); push!(cx3_loc, cz/npts)
        end
        push!(x1mn_loc, x1mn); push!(x1mx_loc, x1mx)
        push!(x2mn_loc, x2mn); push!(x2mx_loc, x2mx)
    end

    gel_g  = MPI.gather(local_gel,  comm)
    cx1_g  = MPI.gather(cx1_loc,    comm)
    cx2_g  = MPI.gather(cx2_loc,    comm)
    cx3_g  = MPI.gather(cx3_loc,    comm)
    x1mn_g = MPI.gather(x1mn_loc,   comm)
    x1mx_g = MPI.gather(x1mx_loc,   comm)
    x2mn_g = MPI.gather(x2mn_loc,   comm)
    x2mx_g = MPI.gather(x2mx_loc,   comm)
    rk_g   = MPI.gather(local_rk,   comm)

    # -----------------------------------------------------------------
    # Phase 2 (rank 0) – build a PAIR LIST: each (child_gidx, parent_gidx,
    # h1, h2, child_is_min_side).  A parent can appear in multiple pairs.
    # -----------------------------------------------------------------
    pair_child_g    = Int[]
    pair_parent_g   = Int[]
    pair_h1_g       = Int[]
    pair_h2_g       = Int[]
    pair_child_ismin = Int[]  # 1 if child is on the min-coordinate side

    if rank == 0
        all_gel  = vcat(gel_g...)
        all_cx1  = vcat(cx1_g...)
        all_cx2  = vcat(cx2_g...)
        all_cx3  = vcat(cx3_g...)
        all_x1mn = vcat(x1mn_g...)
        all_x1mx = vcat(x1mx_g...)
        all_x2mn = vcat(x2mn_g...)
        all_x2mx = vcat(x2mx_g...)
        all_rk   = vcat(rk_g...)

        # Helper: reconstruct (x,y,z) from the (cx1,cx2,cx3) storage convention.
        to_xyz(i) =
            periodic_direction == "periodicx" ? (all_cx3[i], all_cx1[i], all_cx2[i]) :
            periodic_direction == "periodicy" ? (all_cx1[i], all_cx3[i], all_cx2[i]) :
                                               (all_cx1[i], all_cx2[i], all_cx3[i])

        if length(all_gel) > 0
            x3v_min, x3v_max = extrema(all_cx3)
            idx_min = findall(v -> AlmostEqual(v, x3v_min), all_cx3)
            idx_max = findall(v -> AlmostEqual(v, x3v_max), all_cx3)

            # Identify conforming faces by matching centroids on both sides.
            max_cmap = Dict{NTuple{2,Float64},Int}()
            for k in idx_max
                key = (round(all_cx1[k]; digits=5), round(all_cx2[k]; digits=5))
                max_cmap[key] = k
            end
            min_cmap = Dict{NTuple{2,Float64},Int}()
            for k in idx_min
                key = (round(all_cx1[k]; digits=5), round(all_cx2[k]; digits=5))
                min_cmap[key] = k
            end
            conform_keys_min = Set(keys(filter(kv -> haskey(max_cmap, kv[1]), min_cmap)))
            conform_keys_max = Set(keys(filter(kv -> haskey(min_cmap, kv[1]), max_cmap)))

            # Nonconforming parent faces on each side (larger face, no exact centroid match).
            nonconform_max = filter(k -> begin
                ky = (round(all_cx1[k]; digits=5), round(all_cx2[k]; digits=5))
                ky ∉ conform_keys_max
            end, idx_max)
            nonconform_min = filter(k -> begin
                ky = (round(all_cx1[k]; digits=5), round(all_cx2[k]; digits=5))
                ky ∉ conform_keys_min
            end, idx_min)

            tol = 1e-8

            # min-side children (finer on min side, coarser parent on max side).
            # Each min-side nonconforming face is checked against every max-side
            # nonconforming face; up to 4 min children can map to one max parent.
            for k in idx_min
                ky = (round(all_cx1[k]; digits=5), round(all_cx2[k]; digits=5))
                ky in conform_keys_min && continue     # conforming — skip
                c1 = all_cx1[k]; c2 = all_cx2[k]
                for p in nonconform_max
                    c1 > all_x1mn[p]-tol && c1 < all_x1mx[p]+tol &&
                    c2 > all_x2mn[p]-tol && c2 < all_x2mx[p]+tol || continue
                    # k must be the finer (smaller-area) element.
                    (all_x1mx[k]-all_x1mn[k])*(all_x2mx[k]-all_x2mn[k]) <
                    (all_x1mx[p]-all_x1mn[p])*(all_x2mx[p]-all_x2mn[p]) - tol || continue
                    mid1 = (all_x1mn[p]+all_x1mx[p])/2.0
                    mid2 = (all_x2mn[p]+all_x2mx[p])/2.0
                    h1_loc = c1 > mid1 ? 1 : 2
                    h2_loc = c2 > mid2 ? 1 : 2
                    push!(pair_child_g,    k)
                    push!(pair_parent_g,   p)
                    push!(pair_h1_g,       h1_loc)
                    push!(pair_h2_g,       h2_loc)
                    push!(pair_child_ismin, 1)
                    cx,cy,cz = to_xyz(k); px,py,pz = to_xyz(p)
                    # println("[$(periodic_direction)] per-NCF pair (child on min-side): child_gel=$(all_gel[k]) rank=$(all_rk[k]) xyz=($(round(cx;digits=5)),$(round(cy;digits=5)),$(round(cz;digits=5))) | parent_gel=$(all_gel[p]) rank=$(all_rk[p]) xyz=($(round(px;digits=5)),$(round(py;digits=5)),$(round(pz;digits=5))) h1=$(h1_loc) h2=$(h2_loc)")
                    break  # each child has exactly one parent
                end
            end

            # max-side children (finer on max side, coarser parent on min side).
            for k in idx_max
                ky = (round(all_cx1[k]; digits=5), round(all_cx2[k]; digits=5))
                ky in conform_keys_max && continue     # conforming — skip
                c1 = all_cx1[k]; c2 = all_cx2[k]
                for p in nonconform_min
                    c1 > all_x1mn[p]-tol && c1 < all_x1mx[p]+tol &&
                    c2 > all_x2mn[p]-tol && c2 < all_x2mx[p]+tol || continue
                    # k must be the finer (smaller-area) element.
                    (all_x1mx[k]-all_x1mn[k])*(all_x2mx[k]-all_x2mn[k]) <
                    (all_x1mx[p]-all_x1mn[p])*(all_x2mx[p]-all_x2mn[p]) - tol || continue
                    mid1 = (all_x1mn[p]+all_x1mx[p])/2.0
                    mid2 = (all_x2mn[p]+all_x2mx[p])/2.0
                    h1_loc = c1 > mid1 ? 1 : 2
                    h2_loc = c2 > mid2 ? 1 : 2
                    push!(pair_child_g,    k)
                    push!(pair_parent_g,   p)
                    push!(pair_h1_g,       h1_loc)
                    push!(pair_h2_g,       h2_loc)
                    push!(pair_child_ismin, 0)
                    cx,cy,cz = to_xyz(k); px,py,pz = to_xyz(p)
                    # println("[$(periodic_direction)] per-NCF pair (child on max-side): child_gel=$(all_gel[k]) rank=$(all_rk[k]) xyz=($(round(cx;digits=5)),$(round(cy;digits=5)),$(round(cz;digits=5))) | parent_gel=$(all_gel[p]) rank=$(all_rk[p]) xyz=($(round(px;digits=5)),$(round(py;digits=5)),$(round(pz;digits=5))) h1=$(h1_loc) h2=$(h2_loc)")
                    break  # each child has exactly one parent
                end
            end
        end
    end

    # -----------------------------------------------------------------
    # Phase 3 – broadcast the full pair list and the global gel/rank
    # arrays to all ranks so each rank can self-filter.
    # -----------------------------------------------------------------
    MPI.Barrier(comm)

    npairs = (rank == 0) ? Int32(length(pair_child_g)) : Int32(0)
    npairs = MPI.Bcast(npairs, 0, comm)

    total_faces_g = (rank == 0) ? Int32(length(vcat(gel_g...))) : Int32(0)
    total_faces_g = MPI.Bcast(total_faces_g, 0, comm)
    all_gel_bcast = zeros(eltype(elm2pelm), total_faces_g)
    all_rk_bcast  = zeros(Int32, total_faces_g)
    if rank == 0
        all_gel_bcast .= vcat(gel_g...)
        all_rk_bcast  .= vcat(rk_g...)
    end
    MPI.Bcast!(all_gel_bcast, 0, comm)
    MPI.Bcast!(all_rk_bcast,  0, comm)

    pc  = rank == 0 ? pair_child_g    : zeros(Int, npairs)
    pp  = rank == 0 ? pair_parent_g   : zeros(Int, npairs)
    ph1 = rank == 0 ? pair_h1_g       : zeros(Int, npairs)
    ph2 = rank == 0 ? pair_h2_g       : zeros(Int, npairs)
    pm  = rank == 0 ? pair_child_ismin : zeros(Int, npairs)
    MPI.Bcast!(pc,  0, comm)
    MPI.Bcast!(pp,  0, comm)
    MPI.Bcast!(ph1, 0, comm)
    MPI.Bcast!(ph2, 0, comm)
    MPI.Bcast!(pm,  0, comm)

    # Map global element ID → local element index for this rank.
    gel_to_iel = Dict{eltype(elm2pelm), Int}()
    for (i, iel) in enumerate(local_iel)
        gel_to_iel[local_gel[i]] = iel
    end

    # Staging arrays for cross-rank periodic NCF pairs.
    ncf_peri_pg         = Vector{Vector{Int}}()   # [child_iel, child_fid, h1, h2]
    ncf_peri_cg         = Vector{Vector{Int}}()   # [parent_iel, child_fid, h1, h2]
    gpelm_ghost_peri    = eltype(elm2pelm)[]      # parent global element IDs (pg case)
    gpfacets_ghost_peri = Int[]                   # child_fid values (pg case)
    gpfacets_owner_peri = Int[]                   # parent rank (pg case)
    gcelm_ghost_peri    = eltype(elm2pelm)[]      # child global element IDs (cg case)
    gcfacets_ghost_peri = Int[]                   # child_fid values (cg case)
    gcfacets_owner_peri = Int[]                   # child rank (cg case)

    for ipair = 1:npairs
        child_gidx  = pc[ipair];  parent_gidx = pp[ipair]
        h1 = ph1[ipair];          h2 = ph2[ipair]
        ismin = pm[ipair]

        child_gel  = all_gel_bcast[child_gidx]
        parent_gel = all_gel_bcast[parent_gidx]
        R_c = Int(all_rk_bcast[child_gidx])
        R_p = Int(all_rk_bcast[parent_gidx])

        child_fid = ismin == 1 ? fid_min : fid_max

        if R_c == rank && R_p == rank
            # Both local: append directly to non_conforming_facets.
            child_iel  = get(gel_to_iel, child_gel,  -1)
            parent_iel = get(gel_to_iel, parent_gel, -1)
            (child_iel == -1 || parent_iel == -1) && continue
            push!(mesh.non_conforming_facets, [child_iel, parent_iel, child_fid, h1, h2])

        elseif R_c == rank && R_p != rank
            # Child local, parent on remote rank → parent-ghost entry.
            child_iel = get(gel_to_iel, child_gel, -1)
            child_iel == -1 && continue
            push!(ncf_peri_pg,         [child_iel, child_fid, h1, h2])
            push!(gpelm_ghost_peri,    parent_gel)
            push!(gpfacets_ghost_peri, child_fid)
            push!(gpfacets_owner_peri, R_p)

        elseif R_p == rank && R_c != rank
            # Parent local, child on remote rank → child-ghost entry.
            parent_iel = get(gel_to_iel, parent_gel, -1)
            parent_iel == -1 && continue
            push!(ncf_peri_cg,         [parent_iel, child_fid, h1, h2])
            push!(gcelm_ghost_peri,    child_gel)
            push!(gcfacets_ghost_peri, child_fid)
            push!(gcfacets_owner_peri, R_c)
        end
    end

    # Inverse mapping: global element ID → local element index (needed by get_ghost_ips).
    pelm2elm_inv = Dict(elm2pelm[i] => i for i in 1:length(elm2pelm))

    # --- Extend parent-ghost infrastructure for periodic NCF pairs ---
    # Sort by owner rank so get_ghost_ips return order matches our ncf list order.
    if !isempty(gpelm_ghost_peri)
        srt = sortperm(gpfacets_owner_peri)
        gpelm_ghost_peri    = gpelm_ghost_peri[srt]
        gpfacets_ghost_peri = gpfacets_ghost_peri[srt]
        gpfacets_owner_peri = gpfacets_owner_peri[srt]
        ncf_peri_pg         = ncf_peri_pg[srt]
    end
    # All ranks must call get_ghost_ips (MPI collective); empty arrays are fine.
    pgip_peri, pgip_owner_peri = get_ghost_ips(
        gpelm_ghost_peri, gpfacets_ghost_peri, gpfacets_owner_peri,
        mesh.connijk, pelm2elm_inv, mesh.ip2gip, ngl, 1, comm, mesh.SD)
    if !isempty(pgip_peri)
        # println("[DBG pg R$(rank) $(periodic_direction)] get_ghost_ips returned $(length(pgip_peri)) GIPs")
        # println("  pgip_peri     = $(pgip_peri)")
        # println("  pgip_owner    = $(pgip_owner_peri)")
        # println("  gip2ip on THIS rank = $(mesh.gip2ip[pgip_peri])")
        # println("  (0 entries → GIP not owned by this rank)")
    end

    pgip_local_peri = send_and_receive(pgip_peri, pgip_owner_peri, comm)[1]

    mesh.pgip_ghost = vcat(mesh.pgip_ghost, pgip_peri)
    mesh.pgip_owner = vcat(mesh.pgip_owner, pgip_owner_peri)
    mesh.pgip_local = vcat(mesh.pgip_local, pgip_local_peri)
    # Rebuild pgip_local on all ranks (MPI collective).

    for entry in ncf_peri_pg
        push!(mesh.non_conforming_facets_parents_ghost, entry)
        push!(mesh.cip_pg,   entry[1])
        push!(mesh.lfid_pg,  entry[2])
        push!(mesh.half1_pg, entry[3])
        push!(mesh.half2_pg, entry[4])
    end
    # NOTE: do NOT update mesh.num_ncf_pg here — extend_ncf_ip_lists_3D! uses
    # the old interior value to know where to start filling IPc_list_pg.

    # --- Extend child-ghost infrastructure for periodic NCF pairs ---
    if !isempty(gcelm_ghost_peri)
        srt = sortperm(gcfacets_owner_peri)
        gcelm_ghost_peri    = gcelm_ghost_peri[srt]
        gcfacets_ghost_peri = gcfacets_ghost_peri[srt]
        gcfacets_owner_peri = gcfacets_owner_peri[srt]
        ncf_peri_cg         = ncf_peri_cg[srt]
    end
    # All ranks must call get_ghost_ips (MPI collective); empty arrays are fine.
    cgip_peri, cgip_owner_peri = get_ghost_ips(
        gcelm_ghost_peri, gcfacets_ghost_peri, gcfacets_owner_peri,
        mesh.connijk, pelm2elm_inv, mesh.ip2gip, ngl, 2, comm, mesh.SD)
    if !isempty(cgip_peri)
        # println("[DBG cg R$(rank) $(periodic_direction)] get_ghost_ips returned $(length(cgip_peri)) GIPs")
        # println("  cgip_peri     = $(cgip_peri)")
        # println("  cgip_owner    = $(cgip_owner_peri)")
        # println("  gip2ip on THIS rank = $(mesh.gip2ip[cgip_peri])")
        # println("  (0 entries → GIP not owned by this rank)")
    end

    
    cgip_local_peri = send_and_receive( cgip_peri, cgip_owner_peri, comm)[1]

    mesh.cgip_ghost = vcat(mesh.cgip_ghost, cgip_peri)
    mesh.cgip_owner = vcat(mesh.cgip_owner, cgip_owner_peri)
    mesh.cgip_local = vcat(mesh.cgip_local, cgip_local_peri)
    # Rebuild cgip_local on all ranks (MPI collective).

    for entry in ncf_peri_cg
        push!(mesh.non_conforming_facets_children_ghost, entry)
        push!(mesh.pip_cg,   entry[1])
        push!(mesh.lfid_cg,  entry[2])
        push!(mesh.half1_cg, entry[3])
        push!(mesh.half2_cg, entry[4])
    end
    # NOTE: do NOT update mesh.num_ncf_cg here — extend_ncf_ip_lists_3D! uses
    # the old interior value to know where to start filling IPp_list_cg.
end


# ---------------------------------------------------------------------------
# Extend IPc_list, IPp_list (and their ghost variants) after periodic NCF
# pairs have been appended to non_conforming_facets* by
# collect_periodic_ncf_pairs_3D!.  Call once after all three directions.
# ---------------------------------------------------------------------------
function extend_ncf_ip_lists_3D!(mesh)
    ngl  = mesh.ngl
    ngl2 = ngl * ngl

    # Helper: extract ngl² IPs from element iel at face fid.
    function face_ips(iel, fid)
        if fid == 1       # front j=ngl
            return reshape(mesh.connijk[iel, 1:ngl, ngl, 1:ngl], :)
        elseif fid == 2   # back  j=1
            return reshape(mesh.connijk[iel, 1:ngl, 1,   1:ngl], :)
        elseif fid == 3   # bottom k=ngl
            return reshape(mesh.connijk[iel, 1:ngl, 1:ngl, ngl], :)
        elseif fid == 4   # top    k=1
            return reshape(mesh.connijk[iel, 1:ngl, 1:ngl, 1  ], :)
        elseif fid == 5   # right  i=1
            return reshape(mesh.connijk[iel, 1,   1:ngl, 1:ngl], :)
        else              # left   i=ngl
            return reshape(mesh.connijk[iel, ngl, 1:ngl, 1:ngl], :)
        end
    end

    # Opposite face mapping (child's face → parent's face).
    opp = (2, 1, 4, 3, 6, 5)

    # --- Extend IPc_list / IPp_list for local-local pairs ---
    new_num_ncf       = length(mesh.non_conforming_facets)
    old_num_ncf       = mesh.num_ncf
    mesh.num_ncf_peri = new_num_ncf - old_num_ncf
    if new_num_ncf > old_num_ncf
        new_IPc = zeros(Int, ngl2, new_num_ncf)
        new_IPp = zeros(Int, ngl2, new_num_ncf)
        if old_num_ncf > 0
            new_IPc[:, 1:old_num_ncf] .= mesh.IPc_list
            new_IPp[:, 1:old_num_ncf] .= mesh.IPp_list
        end
        for idx = old_num_ncf+1:new_num_ncf
            ciel, piel, cfid, _, _ = mesh.non_conforming_facets[idx]
            new_IPc[:, idx] .= face_ips(ciel, cfid)
            new_IPp[:, idx] .= face_ips(piel, opp[cfid])
        end
        mesh.IPc_list = new_IPc
        mesh.IPp_list = new_IPp
        mesh.num_ncf  = new_num_ncf
    end

    # --- Extend IPc_list_pg for parent-ghost pairs (child local, parent remote) ---
    new_num_pg = length(mesh.non_conforming_facets_parents_ghost)
    old_num_pg = mesh.num_ncf_pg
    mesh.num_ncf_pg_peri = new_num_pg - old_num_pg
    # @info mesh.rank new_num_pg old_num_pg
    if new_num_pg > old_num_pg
        new_IPc_pg = zeros(Int, ngl2, new_num_pg)
        if old_num_pg > 0
            new_IPc_pg[:, 1:old_num_pg] .= mesh.IPc_list_pg
        end
        for idx = old_num_pg+1:new_num_pg
            ciel, cfid, _, _ = mesh.non_conforming_facets_parents_ghost[idx]
            new_IPc_pg[:, idx] .= face_ips(ciel, cfid)
        end
        mesh.IPc_list_pg = new_IPc_pg
        mesh.num_ncf_pg  = new_num_pg
    end

    # --- Extend IPp_list_cg for child-ghost pairs (parent local, child remote) ---
    new_num_cg = length(mesh.non_conforming_facets_children_ghost)
    old_num_cg = mesh.num_ncf_cg
    mesh.num_ncf_cg_peri = new_num_cg - old_num_cg
    if new_num_cg > old_num_cg
        new_IPp_cg = zeros(Int, ngl2, new_num_cg)
        if old_num_cg > 0
            new_IPp_cg[:, 1:old_num_cg] .= mesh.IPp_list_cg
        end
        for idx = old_num_cg+1:new_num_cg
            piel, cfid, _, _ = mesh.non_conforming_facets_children_ghost[idx]
            # Parent is at the opposite face of the child's face (cfid).
            new_IPp_cg[:, idx] .= face_ips(piel, opp[cfid])
        end
        mesh.IPp_list_cg = new_IPp_cg
        mesh.num_ncf_cg  = new_num_cg
    end
end

# ---------------------------------------------------------------------------
# Debug helper: print a detailed table of all periodic NCF entries.
# Call after extend_ncf_ip_lists_3D! so that num_ncf_peri etc. are set.
#
# For each entry prints:
#   rank | type | child iel/gel | child face name | child face centroid |
#   parent iel/gel (or GHOST) | parent face name | h1 h2
#
# Face-ID convention:
#   1 front(j=ngl)  2 back(j=1)  3 bottom(k=ngl)  4 top(k=1)
#   5 right(i=1)    6 left(i=ngl)
# ---------------------------------------------------------------------------
function print_periodic_ncf_debug!(mesh)
    comm  = MPI.COMM_WORLD
    rank  = MPI.Comm_rank(comm)
    ngl   = mesh.ngl

    face_name = ("front(j=ngl)", "back(j=1)", "bottom(k=ngl)", "top(k=1)", "right(i=1)", "left(i=ngl)")
    opp_fid   = (2, 1, 4, 3, 6, 5)

    function direction_from_fid(fid)
        fid in (5, 6) && return "x"
        fid in (1, 2) && return "y"
        return "z"
    end

    function face_centroid(iel, fid)
        if fid == 1
            ips = reshape(mesh.connijk[iel, 1:ngl, ngl, 1:ngl], :)
        elseif fid == 2
            ips = reshape(mesh.connijk[iel, 1:ngl, 1,   1:ngl], :)
        elseif fid == 3
            ips = reshape(mesh.connijk[iel, 1:ngl, 1:ngl, ngl], :)
        elseif fid == 4
            ips = reshape(mesh.connijk[iel, 1:ngl, 1:ngl, 1  ], :)
        elseif fid == 5
            ips = reshape(mesh.connijk[iel, 1,   1:ngl, 1:ngl], :)
        else
            ips = reshape(mesh.connijk[iel, ngl, 1:ngl, 1:ngl], :)
        end
        n  = length(ips)
        cx = sum(mesh.x[ip] for ip in ips) / n
        cy = sum(mesh.y[ip] for ip in ips) / n
        cz = sum(mesh.z[ip] for ip in ips) / n
        return cx, cy, cz
    end

    r5(v) = round(v; digits=5)
    sep = "─────────────────────────────────────────────────────────────────────────"

    # Print a range of local-local NCF entries from non_conforming_facets.
    function print_ll_range(label, i_start, i_end)
        n = i_end - i_start + 1
        println(sep)
        println("[NCF] Rank $rank | $label local-local pairs ($n entries)")
        n == 0 && return
        println("  idx | dir       | child  iel(gel)       | child_face          | child_centroid (x,y,z)")
        println("      |           | parent iel(gel)       | parent_face         | parent_centroid (x,y,z)       | h1 h2")
        for i = i_start:i_end
            ciel, piel, cfid, h1, h2 = mesh.non_conforming_facets[i]
            cgel = mesh.el2gel[ciel]
            pgel = mesh.el2gel[piel]
            pfid = opp_fid[cfid]
            dir  = direction_from_fid(cfid)
            ccx, ccy, ccz = face_centroid(ciel, cfid)
            pcx, pcy, pcz = face_centroid(piel, pfid)
            idx  = i - i_start + 1
            println("  $idx | $(rpad(dir,9)) | child  iel=$(lpad(ciel,5)) gel=$(lpad(cgel,6)) | $(rpad(face_name[cfid],20)) | ($(r5(ccx)), $(r5(ccy)), $(r5(ccz)))")
            println("      |           | parent iel=$(lpad(piel,5)) gel=$(lpad(pgel,6)) | $(rpad(face_name[pfid],20)) | ($(r5(pcx)), $(r5(pcy)), $(r5(pcz))) | h1=$h1 h2=$h2")
        end
    end

    # Print a range of parent-ghost NCF entries (child local, parent remote).
    function print_pg_range(label, i_start, i_end)
        n = i_end - i_start + 1
        println(sep)
        println("[NCF] Rank $rank | $label child-local/parent-ghost pairs ($n entries)")
        n == 0 && return
        println("  idx | dir       | child  iel(gel)       | child_face          | child_centroid (x,y,z)        | h1 h2")
        for i = i_start:i_end
            ciel, cfid, h1, h2 = mesh.non_conforming_facets_parents_ghost[i]
            cgel = mesh.el2gel[ciel]
            pfid = opp_fid[cfid]
            dir  = direction_from_fid(cfid)
            ccx, ccy, ccz = face_centroid(ciel, cfid)
            idx  = i - i_start + 1
            println("  $idx | $(rpad(dir,9)) | child  iel=$(lpad(ciel,5)) gel=$(lpad(cgel,6)) | $(rpad(face_name[cfid],20)) | ($(r5(ccx)), $(r5(ccy)), $(r5(ccz))) | h1=$h1 h2=$h2")
            println("      |           | parent GHOST                    | parent_face=$(face_name[pfid])")
        end
    end

    # Print a range of child-ghost NCF entries (parent local, child remote).
    function print_cg_range(label, i_start, i_end)
        n = i_end - i_start + 1
        println(sep)
        println("[NCF] Rank $rank | $label parent-local/child-ghost pairs ($n entries)")
        n == 0 && return
        println("  idx | dir       | parent iel(gel)       | parent_face         | parent_centroid (x,y,z)       | h1 h2")
        for i = i_start:i_end
            piel, cfid, h1, h2 = mesh.non_conforming_facets_children_ghost[i]
            pgel = mesh.el2gel[piel]
            pfid = opp_fid[cfid]
            dir  = direction_from_fid(cfid)
            pcx, pcy, pcz = face_centroid(piel, pfid)
            idx  = i - i_start + 1
            println("  $idx | $(rpad(dir,9)) | parent iel=$(lpad(piel,5)) gel=$(lpad(pgel,6)) | $(rpad(face_name[pfid],20)) | ($(r5(pcx)), $(r5(pcy)), $(r5(pcz))) | h1=$h1 h2=$h2")
            println("      |           | child  GHOST                    | child_face=$(face_name[cfid])")
        end
    end

    n_tot    = mesh.num_ncf;       n_ll_peri = mesh.num_ncf_peri
    n_pg_tot = mesh.num_ncf_pg;    n_pg_peri = mesh.num_ncf_pg_peri
    n_cg_tot = mesh.num_ncf_cg;    n_cg_peri = mesh.num_ncf_cg_peri

    n_ll_int = n_tot    - n_ll_peri
    n_pg_int = n_pg_tot - n_pg_peri
    n_cg_int = n_cg_tot - n_cg_peri

    print_ll_range("interior ", 1,             n_ll_int)
    print_ll_range("periodic ", n_ll_int + 1,  n_tot)
    print_pg_range("interior ", 1,             n_pg_int)
    print_pg_range("periodic ", n_pg_int + 1,  n_pg_tot)
    print_cg_range("interior ", 1,             n_cg_int)
    print_cg_range("periodic ", n_cg_int + 1,  n_cg_tot)

    println(sep)
    flush(stdout)
end

# ---------------------------------------------------------------------------
# For each NCF pair print the ngl² child IPs and the ngl² parent IPs with
# their (x,y,z) coordinates, indexed so that IPc[j] and IPp[j] share the
# same flat face position j = j1 + (j2-1)*ngl.
#
# Also groups all NCF pairs by parent element and, for each parent face node,
# lists the contributing child nodes across the (up to 4) children — giving
# the "1 parent ip → N child ips" view the DSS gather accumulates.
# ---------------------------------------------------------------------------
function print_ncf_ip_coords!(mesh)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    ngl  = mesh.ngl
    ngl2 = ngl * ngl

    face_name = ("front(j=ngl)", "back(j=1)", "bottom(k=ngl)", "top(k=1)", "right(i=1)", "left(i=ngl)")
    opp_fid   = (2, 1, 4, 3, 6, 5)
    r3(v)     = round(v; digits=3)
    sep = "─────────────────────────────────────────────────────────────────────────"

    n_tot    = mesh.num_ncf
    n_ll_int = n_tot - mesh.num_ncf_peri

    # ── per-pair view: IPc[:] and IPp[:] with coordinates ──────────────────
    println(sep)
    println("[NCF-IPS] Rank $rank | per-pair IPc/IPp coordinates ($n_tot pairs)")
    for idx = 1:n_tot
        ciel, piel, cfid, h1, h2 = mesh.non_conforming_facets[idx]
        pfid  = opp_fid[cfid]
        cgel  = mesh.el2gel[ciel]
        pgel  = mesh.el2gel[piel]
        label = idx <= n_ll_int ? "interior" : "periodic"
        println(sep)
        println("  NCF #$idx ($label)  child iel=$ciel gel=$cgel $(face_name[cfid])  h1=$h1 h2=$h2")
        println("                     parent iel=$piel gel=$pgel $(face_name[pfid])")
        println("    j  | child ip  gip  (x, y, z)                        | parent ip  gip  (x, y, z)")
        for j = 1:ngl2
            cip = mesh.IPc_list[j, idx]
            pip = mesh.IPp_list[j, idx]
            cgip = mesh.ip2gip[cip]; pgip = mesh.ip2gip[pip]
            println("   $(lpad(j,2)) | ip=$(lpad(cip,6)) gip=$(lpad(cgip,6))  ($(r3(mesh.x[cip])),$(r3(mesh.y[cip])),$(r3(mesh.z[cip])))  | ip=$(lpad(pip,6)) gip=$(lpad(pgip,6))  ($(r3(mesh.x[pip])),$(r3(mesh.y[pip])),$(r3(mesh.z[pip])))")
        end
    end

    # ── parent-centric view: for each unique parent element, list all
    #    children and show which child IP maps to each parent IP ─────────────
    println(sep)
    println("[NCF-IPS] Rank $rank | parent-centric view (1 parent ip ← N child ips)")
    # group pair indices by parent iel
    parent_to_pairs = Dict{Int, Vector{Int}}()
    for idx = 1:n_tot
        piel = mesh.non_conforming_facets[idx][2]
        push!(get!(parent_to_pairs, piel, Int[]), idx)
    end
    for (piel, pairs) in sort(collect(parent_to_pairs), by=first)
        pgel = mesh.el2gel[piel]
        println(sep)
        println("  Parent iel=$piel gel=$pgel  ($(length(pairs)) children)")
        # header: parent node | child1 | child2 | ...
        hdr = "  j  | parent ip  gip  (x,y,z)  "
        for idx in pairs
            ciel = mesh.non_conforming_facets[idx][1]
            cfid = mesh.non_conforming_facets[idx][3]
            hdr *= "| child iel=$ciel $(face_name[cfid][1:5])  "
        end
        println(hdr)
        for j = 1:ngl2
            pip  = mesh.IPp_list[j, pairs[1]]
            pgip = mesh.ip2gip[pip]
            row  = "  $(lpad(j,2)) | ip=$(lpad(pip,6)) gip=$(lpad(pgip,6))  ($(r3(mesh.x[pip])),$(r3(mesh.y[pip])),$(r3(mesh.z[pip])))  "
            for idx in pairs
                cip  = mesh.IPc_list[j, idx]
                cgip = mesh.ip2gip[cip]
                row *= "| ip=$(lpad(cip,6)) gip=$(lpad(cgip,6))  ($(r3(mesh.x[cip])),$(r3(mesh.y[cip])),$(r3(mesh.z[cip])))  "
            end
            println(row)
        end
    end
    println(sep)
    flush(stdout)
end

# ---------------------------------------------------------------------------
# Lightweight periodic-NCF *detection* — no NCF interpolation infrastructure.
#
# Runs Phases 1–2 of collect_periodic_ncf_pairs_3D! (centroid-based pair
# detection) and pushes the global element IDs of the coarser-side (parent)
# elements at periodic NCF boundaries into mesh.periodic_ncf_parent_gels.
#
# Called once per periodic direction (periodicx/y/z) inside mod_mesh_read_gmsh!
# when ladaptive == true, replacing the old collect_periodic_ncf_pairs_3D! +
# extend_ncf_ip_lists_3D! calls.  amr_strategy! then iteratively refines the
# detected parents until periodic boundaries are conforming.
# ---------------------------------------------------------------------------
function detect_periodic_ncf_parent_gels!(mesh, periodic_direction, elm2pelm)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    ngl  = mesh.ngl
    ngl2 = ngl * ngl

    # ---- Phase 1: each rank collects its periodic boundary faces ----
    local_iel = Int[]
    local_gel = eltype(elm2pelm)[]
    local_rk  = Int[]
    cx1_loc   = Float64[];  cx2_loc = Float64[];  cx3_loc = Float64[]
    x1mn_loc  = Float64[];  x1mx_loc = Float64[]
    x2mn_loc  = Float64[];  x2mx_loc = Float64[]

    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        mesh.bdy_face_type[iface_bdy] == periodic_direction || continue
        iel = mesh.bdy_face_in_elem[iface_bdy]
        gel = elm2pelm[iel]
        cx = 0.0; cy = 0.0; cz = 0.0
        x1mn = Inf; x1mx = -Inf; x2mn = Inf; x2mx = -Inf
        for kk = 1:ngl, ll = 1:ngl
            ip = mesh.poin_in_bdy_face[iface_bdy, kk, ll]
            xx = mesh.x[ip]; yy = mesh.y[ip]; zz = mesh.z[ip]
            cx += xx; cy += yy; cz += zz
            if periodic_direction == "periodicx"
                x1mn = min(x1mn, yy); x1mx = max(x1mx, yy)
                x2mn = min(x2mn, zz); x2mx = max(x2mx, zz)
            elseif periodic_direction == "periodicy"
                x1mn = min(x1mn, xx); x1mx = max(x1mx, xx)
                x2mn = min(x2mn, zz); x2mx = max(x2mx, zz)
            else
                x1mn = min(x1mn, xx); x1mx = max(x1mx, xx)
                x2mn = min(x2mn, yy); x2mx = max(x2mx, yy)
            end
        end
        npts = Float64(ngl2)
        push!(local_iel, iel); push!(local_gel, gel); push!(local_rk, rank)
        if periodic_direction == "periodicx"
            push!(cx1_loc, cy/npts); push!(cx2_loc, cz/npts); push!(cx3_loc, cx/npts)
        elseif periodic_direction == "periodicy"
            push!(cx1_loc, cx/npts); push!(cx2_loc, cz/npts); push!(cx3_loc, cy/npts)
        else
            push!(cx1_loc, cx/npts); push!(cx2_loc, cy/npts); push!(cx3_loc, cz/npts)
        end
        push!(x1mn_loc, x1mn); push!(x1mx_loc, x1mx)
        push!(x2mn_loc, x2mn); push!(x2mx_loc, x2mx)
    end

    gel_g  = MPI.gather(local_gel,  comm)
    cx1_g  = MPI.gather(cx1_loc,    comm)
    cx2_g  = MPI.gather(cx2_loc,    comm)
    cx3_g  = MPI.gather(cx3_loc,    comm)
    x1mn_g = MPI.gather(x1mn_loc,   comm)
    x1mx_g = MPI.gather(x1mx_loc,   comm)
    x2mn_g = MPI.gather(x2mn_loc,   comm)
    x2mx_g = MPI.gather(x2mx_loc,   comm)
    rk_g   = MPI.gather(local_rk,   comm)

    # ---- Phase 2 (rank 0): detect coarser-side (parent) elements ----
    parent_gels_r0 = eltype(elm2pelm)[]

    if rank == 0
        all_gel  = vcat(gel_g...)
        all_cx1  = vcat(cx1_g...)
        all_cx2  = vcat(cx2_g...)
        all_cx3  = vcat(cx3_g...)
        all_x1mn = vcat(x1mn_g...)
        all_x1mx = vcat(x1mx_g...)
        all_x2mn = vcat(x2mn_g...)
        all_x2mx = vcat(x2mx_g...)

        if length(all_gel) > 0
            x3v_min, x3v_max = extrema(all_cx3)
            idx_min = findall(v -> AlmostEqual(v, x3v_min), all_cx3)
            idx_max = findall(v -> AlmostEqual(v, x3v_max), all_cx3)

            max_cmap = Dict{NTuple{2,Float64},Int}()
            for k in idx_max
                max_cmap[(round(all_cx1[k]; digits=5), round(all_cx2[k]; digits=5))] = k
            end
            min_cmap = Dict{NTuple{2,Float64},Int}()
            for k in idx_min
                min_cmap[(round(all_cx1[k]; digits=5), round(all_cx2[k]; digits=5))] = k
            end
            conform_keys_min = Set(keys(filter(kv -> haskey(max_cmap, kv[1]), min_cmap)))
            conform_keys_max = Set(keys(filter(kv -> haskey(min_cmap, kv[1]), max_cmap)))

            nonconform_max = filter(k -> (round(all_cx1[k]; digits=5), round(all_cx2[k]; digits=5)) ∉ conform_keys_max, idx_max)
            nonconform_min = filter(k -> (round(all_cx1[k]; digits=5), round(all_cx2[k]; digits=5)) ∉ conform_keys_min, idx_min)

            tol = 1e-8
            parent_gel_set = Set{eltype(elm2pelm)}()

            # min-side children → max-side parents
            for k in idx_min
                (round(all_cx1[k]; digits=5), round(all_cx2[k]; digits=5)) in conform_keys_min && continue
                c1 = all_cx1[k]; c2 = all_cx2[k]
                for p in nonconform_max
                    c1 > all_x1mn[p]-tol && c1 < all_x1mx[p]+tol &&
                    c2 > all_x2mn[p]-tol && c2 < all_x2mx[p]+tol || continue
                    (all_x1mx[k]-all_x1mn[k])*(all_x2mx[k]-all_x2mn[k]) <
                    (all_x1mx[p]-all_x1mn[p])*(all_x2mx[p]-all_x2mn[p]) - tol || continue
                    push!(parent_gel_set, all_gel[p])
                    break
                end
            end
            # max-side children → min-side parents
            for k in idx_max
                (round(all_cx1[k]; digits=5), round(all_cx2[k]; digits=5)) in conform_keys_max && continue
                c1 = all_cx1[k]; c2 = all_cx2[k]
                for p in nonconform_min
                    c1 > all_x1mn[p]-tol && c1 < all_x1mx[p]+tol &&
                    c2 > all_x2mn[p]-tol && c2 < all_x2mx[p]+tol || continue
                    (all_x1mx[k]-all_x1mn[k])*(all_x2mx[k]-all_x2mn[k]) <
                    (all_x1mx[p]-all_x1mn[p])*(all_x2mx[p]-all_x2mn[p]) - tol || continue
                    push!(parent_gel_set, all_gel[p])
                    break
                end
            end
            parent_gels_r0 = collect(parent_gel_set)
        end
    end

    # ---- Phase 3: broadcast parent gel IDs; each rank records locally-owned ones ----
    nparents = (rank == 0) ? Int32(length(parent_gels_r0)) : Int32(0)
    nparents = MPI.Bcast(nparents, 0, comm)
    parent_gels_all = zeros(eltype(elm2pelm), nparents)
    if rank == 0
        parent_gels_all .= parent_gels_r0
    end
    MPI.Bcast!(parent_gels_all, 0, comm)

    gel_to_iel = Dict{eltype(elm2pelm), Int}(local_gel[i] => local_iel[i] for i in eachindex(local_iel))
    for gel in parent_gels_all
        if haskey(gel_to_iel, gel)
            push!(mesh.periodic_ncf_parent_gels, Int64(gel))
        end
    end
end

function detect_periodic_ncf_parent_gels_2D!(mesh, periodic_direction, elm2pelm)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    ngl  = mesh.ngl

    # ---- Phase 1: each rank collects its periodic boundary edges ----
    local_iel = Int[]
    local_gel = eltype(elm2pelm)[]
    local_rk  = Int[]
    cx1_loc   = Float64[]   # tangential centroid coordinate (along the edge)
    cx3_loc   = Float64[]   # normal centroid coordinate (discriminates the two sides)
    x1mn_loc  = Float64[]   # edge extent (tangential min)
    x1mx_loc  = Float64[]   # edge extent (tangential max)

    for iedge_bdy = 1:size(mesh.bdy_edge_type, 1)
        mesh.bdy_edge_type[iedge_bdy] == periodic_direction || continue
        iel = mesh.bdy_edge_in_elem[iedge_bdy]
        gel = elm2pelm[iel]
        cx = 0.0; cy = 0.0
        x1mn = Inf; x1mx = -Inf
        for k = 1:ngl
            ip = mesh.poin_in_bdy_edge[iedge_bdy, k]
            xx = mesh.x[ip]; yy = mesh.y[ip]
            cx += xx; cy += yy
            if periodic_direction == "periodicx"
                x1mn = min(x1mn, yy); x1mx = max(x1mx, yy)
            else  # "periodicy"
                x1mn = min(x1mn, xx); x1mx = max(x1mx, xx)
            end
        end
        npts = Float64(ngl)
        push!(local_iel, iel); push!(local_gel, gel); push!(local_rk, rank)
        if periodic_direction == "periodicx"
            push!(cx1_loc, cy/npts); push!(cx3_loc, cx/npts)
        else  # "periodicy"
            push!(cx1_loc, cx/npts); push!(cx3_loc, cy/npts)
        end
        push!(x1mn_loc, x1mn); push!(x1mx_loc, x1mx)
    end

    gel_g  = MPI.gather(local_gel, comm)
    cx1_g  = MPI.gather(cx1_loc,  comm)
    cx3_g  = MPI.gather(cx3_loc,  comm)
    x1mn_g = MPI.gather(x1mn_loc, comm)
    x1mx_g = MPI.gather(x1mx_loc, comm)

    # ---- Phase 2 (rank 0): detect coarser-side (parent) elements ----
    parent_gels_r0 = eltype(elm2pelm)[]

    if rank == 0
        all_gel  = vcat(gel_g...)
        all_cx1  = vcat(cx1_g...)
        all_cx3  = vcat(cx3_g...)
        all_x1mn = vcat(x1mn_g...)
        all_x1mx = vcat(x1mx_g...)

        if length(all_gel) > 0
            x3v_min, x3v_max = extrema(all_cx3)
            idx_min = findall(v -> AlmostEqual(v, x3v_min), all_cx3)
            idx_max = findall(v -> AlmostEqual(v, x3v_max), all_cx3)

            # Use sigdigits (relative precision) so that small centroid values near
            # z=0 in stretched meshes are still distinguished correctly.
            # digits=5 (absolute) would map both 1.49e-5 and 7.45e-6 to 1e-5,
            # causing false conforming matches in the NCF detection.
            key(v) = round(v; sigdigits=6)

            max_cmap = Dict{Float64, Int}()
            for k in idx_max
                max_cmap[key(all_cx1[k])] = k
            end
            min_cmap = Dict{Float64, Int}()
            for k in idx_min
                min_cmap[key(all_cx1[k])] = k
            end
            conform_keys_min = Set(keys(filter(kv -> haskey(max_cmap, kv[1]), min_cmap)))
            conform_keys_max = Set(keys(filter(kv -> haskey(min_cmap, kv[1]), max_cmap)))

            nonconform_max = filter(k -> key(all_cx1[k]) ∉ conform_keys_max, idx_max)
            nonconform_min = filter(k -> key(all_cx1[k]) ∉ conform_keys_min, idx_min)

            tol = 1e-8
            parent_gel_set = Set{eltype(elm2pelm)}()

            # min-side children → max-side parents
            for k in idx_min
                key(all_cx1[k]) in conform_keys_min && continue
                c1 = all_cx1[k]
                for p in nonconform_max
                    c1 > all_x1mn[p]-tol && c1 < all_x1mx[p]+tol || continue
                    (all_x1mx[k]-all_x1mn[k]) < (all_x1mx[p]-all_x1mn[p]) - tol || continue
                    push!(parent_gel_set, all_gel[p])
                    break
                end
            end
            # max-side children → min-side parents
            for k in idx_max
                key(all_cx1[k]) in conform_keys_max && continue
                c1 = all_cx1[k]
                for p in nonconform_min
                    c1 > all_x1mn[p]-tol && c1 < all_x1mx[p]+tol || continue
                    (all_x1mx[k]-all_x1mn[k]) < (all_x1mx[p]-all_x1mn[p]) - tol || continue
                    push!(parent_gel_set, all_gel[p])
                    break
                end
            end
            parent_gels_r0 = collect(parent_gel_set)
        end
    end

    # ---- Phase 3: broadcast parent gel IDs; each rank records locally-owned ones ----
    nparents = (rank == 0) ? Int32(length(parent_gels_r0)) : Int32(0)
    nparents = MPI.Bcast(nparents, 0, comm)
    parent_gels_all = zeros(eltype(elm2pelm), nparents)
    if rank == 0
        parent_gels_all .= parent_gels_r0
    end
    MPI.Bcast!(parent_gels_all, 0, comm)

    gel_to_iel = Dict{eltype(elm2pelm), Int}(local_gel[i] => local_iel[i] for i in eachindex(local_iel))
    for gel in parent_gels_all
        if haskey(gel_to_iel, gel)
            push!(mesh.periodic_ncf_parent_gels, Int64(gel))
        end
    end
end
