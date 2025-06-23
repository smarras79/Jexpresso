function periodicity_restructure!(mesh,
                                  coords,
                                  xmax,xmin,ymax,ymin,zmax,zmin,poin_in_bdy_face,poin_in_bdy_edge,ngl,ngr,nelem,npoin,nsd,bdy_edge_type,
        bdy_face_type,bdy_face_in_elem,bdy_edge_in_elem,connijk,connijk_lag,npoin_linear,nelem_semi_inf,
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
