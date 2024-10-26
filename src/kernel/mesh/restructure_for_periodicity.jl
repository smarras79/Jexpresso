function periodicity_restructure!(mesh,inputs,backend)
    #per1 = inputs[:per1dd]
    #per2 = inputs[:per2]
    #determine boundary vectors
    if (mesh.nsd == 2)
	if (inputs[:lperiodic_laguerre] && "Laguerre" in mesh.bdy_edge_type)
            xmin = minimum(mesh.x)
            xmax = maximum(mesh.x)
            e1 = 0
            e = 1
            i1 = 0
            i2 = 0
            while (e1 == 0)
                ip_temp = mesh.connijk_lag[e,1,1]
                ip_temp1 = mesh.connijk_lag[e,mesh.ngl,1]
                if (AlmostEqual(mesh.x[ip_temp],xmin))
                    e1 = e
                    i1 = 1
                end
                if (AlmostEqual(mesh.x[ip_temp1],xmin))
                    e1 = e
                    i1 = mesh.ngl
                end
                e += 1

            end
            e2 = 0
            e = 1
            while (e2 == 0)
                ip_temp = mesh.connijk_lag[e,1,1]
                ip_temp1 = mesh.connijk_lag[e,mesh.ngl,1]
                if (AlmostEqual(mesh.x[ip_temp],xmax))
                    e2 = e
                    i2 = 1
                end
                if (AlmostEqual(mesh.x[ip_temp1],xmax))
                    e2 = e
                    i2 = mesh.ngl
                end
                e += 1

            end
            for j = 1:mesh.ngr
                ip1 = mesh.connijk_lag[e1,i1,j]
                ip2 = mesh.connijk_lag[e2,i2,j]
                if (mesh.x[ip1] < mesh.x[ip2])
                    ip_dest = ip1
                    ip_kill = ip2
                elseif(mesh.x[ip1] > mesh.x[ip2])
                    ip_dest = ip2
                    ip_kill = ip1
                end
                mesh.connijk_lag[e1,i1,j] = ip_dest
                mesh.connijk_lag[e2,i2,j] = ip_dest
	        if (j > 1)
                    for ip = ip_kill:mesh.npoin-1
                        mesh.x[ip] = mesh.x[ip+1]
                        mesh.y[ip] = mesh.y[ip+1]
                    end
                    for e=1:mesh.nelem_semi_inf
                        for i=1:mesh.ngl
                            for j=1:mesh.ngr
                                ip = mesh.connijk_lag[e,i,j]
                                if (ip > ip_kill)
                                    mesh.connijk_lag[e,i,j] = ip-1
                                end
                            end
                        end
                    end
                else
                    #  mesh.x[ip_dest] = xmin
                end
            end
            mesh.npoin -= (mesh.ngr-1)
        end  
        if ("periodic1" in mesh.bdy_edge_type)
            finder = false
            iedge_bdy = 1
            while (finder == false)
                if (mesh.bdy_edge_type[iedge_bdy] == "periodic1")
                    ip = mesh.poin_in_bdy_edge[iedge_bdy,1]
                    ip1 = mesh.poin_in_bdy_edge[iedge_bdy,2]
                    per1 = [mesh.x[ip] - mesh.x[ip1],mesh.y[ip] - mesh.y[ip1]]
                    finder = true
                else
                    iedge_bdy +=1
                end
            end
        else
            per1 = [0.0, 1.0]
        end
        if ("periodic2" in mesh.bdy_edge_type)
            finder = false
            iedge_bdy = 1
            while (finder == false)
                if (mesh.bdy_edge_type[iedge_bdy] == "periodic2")
                    ip = mesh.poin_in_bdy_edge[iedge_bdy,1]
                    ip1 = mesh.poin_in_bdy_edge[iedge_bdy,2]
                    per2 = [mesh.x[ip] - mesh.x[ip1],mesh.y[ip] - mesh.y[ip1]]
                    finder = true
                else
                    iedge_bdy +=1
                end
            end 
        else
            per2 = [1.0, 0.0]
        end     
        xx = zeros(TFloat, size(mesh.x,1),1)#KernelAbstractions.zeros(CPU(), TFloat, size(mesh.x,1),1)
        yy = zeros(TFloat, size(mesh.y,1),1)#KernelAbstractions.zeros(CPU(), TFloat, size(mesh.y,1),1)
        poin_bdy= zeros(TInt, size(mesh.poin_in_bdy_edge))#KernelAbstractions.zeros(CPU(), TInt,size(mesh.poin_in_bdy_edge))
        xx .= mesh.x
        yy .= mesh.y
        poin_bdy .=mesh.poin_in_bdy_edge
        interval = [2,3,4]
        for iedge_bdy =1:size(mesh.bdy_edge_type,1)
            if ("periodic2" in mesh.bdy_edge_type && "periodic1" in mesh.bdy_edge_type)
                iel = mesh.bdy_edge_in_elem[iedge_bdy]
                for k=1:mesh.ngl
                    ip = mesh.poin_in_bdy_edge[iedge_bdy,k]
                    if (ip in interval)
                        ip_kill=ip
                        ip_dest=1
                        for e=1:mesh.nelem
                            for ii=1:mesh.ngl
                                for jj=1:mesh.ngl
                                    if (mesh.connijk[iel,ii,jj] == ip_kill)
                                        mesh.connijk[iel,ii,jj] = ip_dest
                                    end
                                end
                            end
                        end
                        for i=ip_kill:mesh.npoin-1
                            mesh.x[i] = mesh.x[i+1]
                            mesh.y[i] = mesh.y[i+1]
                        end
                        mesh.npoin = mesh.npoin-1
                        for iedge =1:size(mesh.poin_in_bdy_edge,1)
                            for kk=1:mesh.ngl
                                val = mesh.poin_in_bdy_edge[iedge,kk]
                                if (val > ip_kill)
                                    mesh.poin_in_bdy_edge[iedge,kk] = val - 1
                                end
                                if (val == ip_kill)
                                    mesh.poin_in_bdy_edge[iedge,kk] = ip_dest
                                end
                            end
                        end

                        for e=1:mesh.nelem
                            for ii=1:mesh.ngl
                                for jj=1:mesh.ngl
                                    ipp = mesh.connijk[e,ii,jj]
                                    if (ipp > ip_kill)
                                        mesh.connijk[e,ii,jj] -= 1
                                    end
                                end
                            end
                        end
                        if (inputs[:lperiodic_laguerre])
                            for e=1:mesh.nelem_semi_inf
                                for ii=1:mesh.ngl
                                    for jj=1:mesh.ngr
                                        ipp = mesh.connijk_lag[e,ii,jj]
                                        if (ipp > ip_kill)
                                            mesh.connijk_lag[e,ii,jj] -= 1
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
        for iedge_bdy =1:size(mesh.bdy_edge_type,1)
            if (mesh.bdy_edge_type[iedge_bdy] == "periodic1" || mesh.bdy_edge_type[iedge_bdy] == "periodic2")
                iel = mesh.bdy_edge_in_elem[iedge_bdy]
                for k =1:mesh.ngl
                    ip = poin_bdy[iedge_bdy, k]
                    ip_true = mesh.poin_in_bdy_edge[iedge_bdy,k]
                    x1 = xx[ip]
                    y1 = yy[ip]
                    m=1
                    l=1
                    for ii=1:mesh.ngl
                        for jj=1:mesh.ngl
                            if (mesh.connijk[iel,ii,jj] == ip_true)
                                l=ii
                                m=jj
                            end
                        end
                    end
                    if (k < mesh.ngl)
                        ip1 = poin_bdy[iedge_bdy,k+1]
                    else
                        ip1 = poin_bdy[iedge_bdy,k-1]
                    end
                    x3 = xx[ip1]
                    y3 = yy[ip1]
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
                    for iedge_per = iedge_bdy+1:size(mesh.bdy_edge_type,1)
                        if (mesh.bdy_edge_type[iedge_per] == mesh.bdy_edge_type[iedge_bdy])
                            iel_per = mesh.bdy_edge_in_elem[iedge_per]
                            for k_per=1:mesh.ngl
                                ip_per = poin_bdy[iedge_per,k_per]
                                ip_true1 = mesh.poin_in_bdy_edge[iedge_per,k_per]
                                x2 = xx[ip_per]
                                y2 = yy[ip_per]
                                if (k_per < mesh.ngl)
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
                                    for ii=1:mesh.ngl
                                        for jj=1:mesh.ngl
                                            if (mesh.connijk[iel_per,ii,jj] == ip_true1)
                                                l1=ii
                                                m1=jj
                                            end
                                        end
                                    end

                                    if (mesh.x[ip_true] < mesh.x[ip_true1])
                                        ip_dest = ip_true
                                        ip_kill = ip_true1
                                    elseif (mesh.x[ip_true] >= mesh.x[ip_true1])
                                        ip_dest = ip_true1
                                        ip_kill = ip_true
                                    end
                                    if (mesh.x[ip_true]*abs(mesh.y[ip_true]) < mesh.x[ip_true1]*abs(mesh.y[ip_true1]) ||mesh.y[ip_true]*abs(mesh.x[ip_true]) < mesh.y[ip_true1]*abs(mesh.x[ip_true1]) )
                                        ip_dest = ip_true
                                        ip_kill = ip_true1
                                    else
                                        ip_dest = ip_true1
                                        ip_kill = ip_true
                                    end 
                                    #@info ip_kill, ip_dest, mesh.x[ip_kill],mesh.x[ip_dest],mesh.y[ip_kill],mesh.y[ip_dest]
                                    mesh.connijk[iel_per,l1,m1] = ip_dest
                                    mesh.connijk[iel,l,m] = ip_dest
                                    mesh.poin_in_bdy_edge[iedge_per,k_per] = ip_dest
                                    mesh.poin_in_bdy_edge[iedge_bdy,k] = ip_dest
                                    ip_true = ip_dest
                                    if !(ip_kill in mesh.connijk)

                                        for i=ip_kill:mesh.npoin-1
                                            mesh.x[i] = mesh.x[i+1]
                                            mesh.y[i] = mesh.y[i+1]
                                        end
                                        mesh.npoin = mesh.npoin-1
                                        for iedge =1:size(mesh.poin_in_bdy_edge,1)
                                            for kk=1:mesh.ngl
                                                val = mesh.poin_in_bdy_edge[iedge,kk]
                                                if (val > ip_kill)
                                                    mesh.poin_in_bdy_edge[iedge,kk] = val - 1
                                                end
                                                if (val == ip_kill)
                                                    mesh.poin_in_bdy_edge[iedge,kk] = ip_dest
                                                end
                                            end
                                        end

                                        for e=1:mesh.nelem
                                            for ii=1:mesh.ngl
                                                for jj=1:mesh.ngl
                                                    ipp = mesh.connijk[e,ii,jj]
                                                    if (ipp > ip_kill)
                                                        mesh.connijk[e,ii,jj] -= 1
                                                    end
                                                    if (ipp == ip_kill)
                                                        mesh.connijk[e,ii,jj] = ip_dest
                                                    end
                                                end
                                            end
                                        end
                                        if (inputs[:lperiodic_laguerre])
                                            for e=1:mesh.nelem_semi_inf
                                                for ii=1:mesh.ngl
                                                    for jj=1:mesh.ngr
                                                        ipp = mesh.connijk_lag[e,ii,jj]
                                                        if (ipp > ip_kill)
                                                            mesh.connijk_lag[e,ii,jj] -= 1
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
    elseif (mesh.nsd == 3)
        if ("periodic1" in mesh.bdy_face_type)
            finder = false
            iface_bdy = 1
            while (finder == false)
                if (mesh.bdy_face_type[iface_bdy] == "periodic1")
                    ip = mesh.poin_in_bdy_face[iface_bdy,1,1]
                    ip1 = mesh.poin_in_bdy_face[iface_bdy,1,2]
                    ip2 = mesh.poin_in_bdy_face[iface_bdy,2,1]
                    t1 = [mesh.x[ip] - mesh.x[ip1],mesh.y[ip] - mesh.y[ip1], mesh.z[ip] - mesh.z[ip1]]
                    t2 = [mesh.x[ip] - mesh.x[ip2],mesh.y[ip] - mesh.y[ip2], mesh.z[ip] - mesh.z[ip2]]
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
        if ("periodic2" in mesh.bdy_face_type)
            finder = false
            iface_bdy = 1
            while (finder == false)
                if (mesh.bdy_face_type[iface_bdy] == "periodic2")
                    ip = mesh.poin_in_bdy_face[iface_bdy,1,1]
                    ip1 = mesh.poin_in_bdy_face[iface_bdy,1,2]
                    ip2 = mesh.poin_in_bdy_face[iface_bdy,2,1]
                    t1 = [mesh.x[ip] - mesh.x[ip1],mesh.y[ip] - mesh.y[ip1], mesh.z[ip] - mesh.z[ip1]]
                    t2 = [mesh.x[ip] - mesh.x[ip2],mesh.y[ip] - mesh.y[ip2], mesh.z[ip] - mesh.z[ip2]]
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
        if ("periodic3" in mesh.bdy_face_type)
            finder = false
            iface_bdy = 1
            while (finder == false)
                if (mesh.bdy_face_type[iface_bdy] == "periodic3")
                    ip = mesh.poin_in_bdy_face[iface_bdy,1,1]
                    ip1 = mesh.poin_in_bdy_face[iface_bdy,1,2]
                    ip2 = mesh.poin_in_bdy_face[iface_bdy,2,1]
                    t1 = [mesh.x[ip] - mesh.x[ip1],mesh.y[ip] - mesh.y[ip1], mesh.z[ip] - mesh.z[ip1]]
                    t2 = [mesh.x[ip] - mesh.x[ip2],mesh.y[ip] - mesh.y[ip2], mesh.z[ip] - mesh.z[ip2]]
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

        xx = zeros(TFloat, size(mesh.x,1),1)#KernelAbstractions.zeros(CPU(), TFloat, size(mesh.x,1),1)
        yy = zeros(TFloat, size(mesh.y,1),1)#KernelAbstractions.zeros(CPU(), TFloat, size(mesh.y,1),1)
        zz = zeros(TFloat, size(mesh.y,1),1)
        poin_bdy= zeros(TInt, size(mesh.poin_in_bdy_face))#KernelAbstractions.zeros(CPU(), TInt,size(mesh.poin_in_bdy_edge))
        xx .= mesh.x
        yy .= mesh.y
        zz .= mesh.z
        poin_bdy .=mesh.poin_in_bdy_face
        ### triple periodicity for corners 
        interval = [2,3,4,5,6,7,8]
        if ("periodic1" in mesh.bdy_face_type && "periodic2" in mesh.bdy_face_type && "periodic3" in mesh.bdy_face_type)    
            for iface_bdy =1:size(mesh.bdy_face_type,1)
                iel = mesh.bdy_face_in_elem[iface_bdy]
                for k=1:mesh.ngl
                    for l=1:mesh.ngl
                        ip = mesh.poin_in_bdy_face[iface_bdy,k,l]
                        if (ip in interval)
                            ip_kill=ip
                            ip_dest=1
                            for e=1:mesh.nelem
                                for ii=1:mesh.ngl
                                    for jj=1:mesh.ngl
                                        for kk=1:mesh.ngl
                                            if (mesh.connijk[iel,ii,jj,kk] == ip_kill)
                                                mesh.connijk[iel,ii,jj,kk] = ip_dest
                                            end
                                        end
                                    end
                                end
                            end
                            for i=ip_kill:mesh.npoin-1
                                mesh.x[i] = mesh.x[i+1]
                                mesh.y[i] = mesh.y[i+1]
                                mesh.z[i] = mesh.z[i+1]
                            end
                            mesh.npoin = mesh.npoin-1
                            for iface =1:size(mesh.poin_in_bdy_face,1)
                                for kk=1:mesh.ngl
                                    for ll=1:mesh.ngl
                                        val = mesh.poin_in_bdy_face[iface,kk,ll]
                                        if (val > ip_kill)
                                            mesh.poin_in_bdy_face[iface,kk,ll] = val - 1
                                        end
                                        if (val == ip_kill)
                                            mesh.poin_in_bdy_face[iface,kk,ll] = ip_dest
                                        end
                                    end
                                end
                            end

                            for e=1:mesh.nelem
                                for ii=1:mesh.ngl
                                    for jj=1:mesh.ngl
                                        for kk=1:mesh.ngl
                                            ipp = mesh.connijk[e,ii,jj,kk]
                                            if (ipp > ip_kill)
                                                mesh.connijk[e,ii,jj,kk] -= 1
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
        if ("periodic1" in mesh.bdy_face_type)
            nperiodic +=1
            double1 = nor1*(mesh.xmax - mesh.xmin)
        end
        if ("periodic2" in mesh.bdy_face_type)
            nperiodic +=1
            double2 = nor2*(mesh.zmax - mesh.zmin)
        end
        if ("periodic3" in mesh.bdy_face_type)
            nperiodic +=1
            if (double1 == [0 ,0 ,0])
                double1 = nor3*(mesh.ymax-mesh.ymin)
            end
            if (double2 == [0 ,0 ,0])
                double2 = nor3*(mesh.ymax-mesh.ymin)
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
                vec = [mesh.x[1] - mesh.x[i], mesh.y[1] - mesh.y[i], mesh.z[1] - mesh.z[i]]
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
            for i in interval
                ip_kill=i
                ip_dest=target_idx
                for e=1:mesh.nelem
                    for ii=1:mesh.ngl
                        for jj=1:mesh.ngl
                            for kk=1:mesh.ngl
                                if (mesh.connijk[e,ii,jj,kk] == ip_kill)
                                    mesh.connijk[e,ii,jj,kk] = ip_dest
                                end
                            end
                        end
                    end
                end
                for i=ip_kill:mesh.npoin-1
                    mesh.x[i] = mesh.x[i+1]
                    mesh.y[i] = mesh.y[i+1]
                    mesh.z[i] = mesh.z[i+1]
                end
                mesh.npoin = mesh.npoin-1
                for iface =1:size(mesh.poin_in_bdy_face,1)
                    for kk=1:mesh.ngl
                        for ll=1:mesh.ngl
                            val = mesh.poin_in_bdy_face[iface,kk,ll]
                            if (val > ip_kill)
                                mesh.poin_in_bdy_face[iface,kk,ll] = val - 1
                            end
                            if (val == ip_kill)
                                mesh.poin_in_bdy_face[iface,kk,ll] = ip_dest
                            end
                        end
                    end
                end
                for e=1:mesh.nelem
                    for ii=1:mesh.ngl
                        for jj=1:mesh.ngl
                            for kk=1:mesh.ngl
                                ipp = mesh.connijk[e,ii,jj,kk]
                                if (ipp > ip_kill)
                                    mesh.connijk[e,ii,jj,kk] -= 1
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
            #make sure to pick a corner consistent with the vtk unwrap
            for i=2:size(plane2,1)
                ii = plane2[i]
                xt = mesh.x[target_idx]
                yt = mesh.y[target_idx]
                zt = mesh.z[target_idx]
                xi = mesh.x[ii]
                yi = mesh.y[ii]
                zi = mesh.z[ii]
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
                #@info xi,yi,zi,xt,yt,zt, comp1,comp2,comp3
                if (comp1 || comp2 || comp3)#(comp1 && comp2) || (comp1 && comp3) || (comp2 && comp3)
                    target_idx = plane2[i]
                    p2_idx = i
                end
            end
            @info mesh.x[target_idx], mesh.y[target_idx],mesh.z[target_idx]
            if !(plane2[1] == target_idx)
                aux = plane2[1]
                plane2[1] = target_idx
                plane2[p2_idx] = aux
            end
            interval = [0, 0, 0]
            for i =1:3
                interval[i] = plane2[i+1]
            end
            for i in interval
                ip_kill=i
                ip_dest=target_idx
                for e=1:mesh.nelem
                    for ii=1:mesh.ngl
                        for jj=1:mesh.ngl
                            for kk=1:mesh.ngl
                                if (mesh.connijk[e,ii,jj,kk] == ip_kill)
                                    mesh.connijk[e,ii,jj,kk] = ip_dest
                                end
                            end
                        end
                    end
                end
                for i=ip_kill:mesh.npoin-1
                    mesh.x[i] = mesh.x[i+1]
                    mesh.y[i] = mesh.y[i+1]
                    mesh.z[i] = mesh.z[i+1]
                end
                mesh.npoin = mesh.npoin-1
                for iface =1:size(mesh.poin_in_bdy_face,1)
                    for kk=1:mesh.ngl
                        for ll=1:mesh.ngl
                            val = mesh.poin_in_bdy_face[iface,kk,ll]
                            if (val > ip_kill)
                                mesh.poin_in_bdy_face[iface,kk,ll] = val - 1
                            end
                            if (val == ip_kill)
                                mesh.poin_in_bdy_face[iface,kk,ll] = ip_dest
                            end
                        end
                    end
                end
                for e=1:mesh.nelem
                    for ii=1:mesh.ngl
                        for jj=1:mesh.ngl
                            for kk=1:mesh.ngl
                                ipp = mesh.connijk[e,ii,jj,kk]
                                if (ipp > ip_kill)
                                    mesh.connijk[e,ii,jj,kk] -= 1
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
        for iface_bdy =1:size(mesh.bdy_face_type,1)
            for k=1:mesh.ngl
                for l=1:mesh.ngl
                    ip = mesh.poin_in_bdy_face[iface_bdy,k,l]
                    if (mesh.bdy_face_type[iface_bdy] == "periodic1")
                        per1_points = [per1_points; ip]
                    elseif (mesh.bdy_face_type[iface_bdy] == "periodic2")
                        per2_points = [per2_points; ip]
                    elseif (mesh.bdy_face_type[iface_bdy] == "periodic3")
                        per3_points = [per3_points; ip]
                    end
                end
            end
        end
        ### remove duplicates
        unique!(per1_points)
        unique!(per2_points)
        unique!(per3_points)
        ### work single periodicity on the periodic1 boundaries if applicable
        if (size(per1_points,1) > 1) 
            for i=1:size(per1_points,1)
                found = false
                i1 = i+1
                while (i1 <= size(per1_points,1) && found == false)
                    ip = per1_points[i]
                    ip1 = per1_points[i1]
                    vec = [mesh.x[ip] - mesh.x[ip1], mesh.y[ip] - mesh.y[ip1], mesh.z[ip] - mesh.z[ip1]]
                    if (determine_colinearity(vec, nor1))
                        #@info "found a match", vec, nor1, ip, ip1, mesh.x[ip], mesh.x[ip1]
                        found = true
                        cond1 = mesh.x[ip]*abs(mesh.y[ip])*abs(mesh.z[ip]) < mesh.x[ip1]*abs(mesh.y[ip1])*abs(mesh.z[ip1])
                        cond2 = mesh.y[ip]*abs(mesh.x[ip])*abs(mesh.z[ip]) < mesh.y[ip1]*abs(mesh.x[ip1])*abs(mesh.z[ip1])
                        cond3 = mesh.z[ip]*abs(mesh.x[ip])*abs(mesh.y[ip]) < mesh.z[ip1]*abs(mesh.x[ip1])*abs(mesh.y[ip1])
                        if (cond1 || cond2 || cond3)    
                            ip_dest = ip
                            ip_kill = ip1
                        else
                            ip_dest = ip1
                            ip_kill = ip
                        end
                        for e=1:mesh.nelem
                            for ii=1:mesh.ngl
                                for jj=1:mesh.ngl
                                    for kk=1:mesh.ngl
                                        ipp = mesh.connijk[e,ii,jj,kk]
                                        if (ipp == ip1 || ipp == ip)
                                            mesh.connijk[e,ii,jj,kk] = ip_dest
                                        end
                                    end
                                end
                            end
                        end
                        for iface = 1:size(mesh.poin_in_bdy_face,1)
                            for kk =1:mesh.ngl
                                for ll=1:mesh.ngl
                                    ipp = mesh.poin_in_bdy_face[iface,kk,ll]
                                    if (ipp == ip1 || ipp == ip)
                                        mesh.poin_in_bdy_face[iface,kk,ll] = ip_dest
                                    end
                                end
                            end
                        end
                        if !(ip_kill in mesh.connijk)
                            mesh.x[ip_dest] = min(mesh.x[ip_kill],mesh.x[ip_dest])
                            mesh.y[ip_dest] = min(mesh.y[ip_kill],mesh.y[ip_dest])
                            mesh.z[ip_dest] = min(mesh.z[ip_kill],mesh.z[ip_dest])
                            for ii=ip_kill:mesh.npoin-1
                                mesh.x[ii] = mesh.x[ii+1]
                                mesh.y[ii] = mesh.y[ii+1]
                                mesh.z[ii] = mesh.z[ii+1]
                            end
                            mesh.npoin = mesh.npoin-1
                            for iface =1:size(mesh.poin_in_bdy_face,1)
                                for kk=1:mesh.ngl
                                    for ll=1:mesh.ngl
                                        val = mesh.poin_in_bdy_face[iface,kk,ll]
                                        if (val > ip_kill)
                                            mesh.poin_in_bdy_face[iface,kk,ll] = val - 1
                                        end
                                        if (val == ip_kill)
                                            mesh.poin_in_bdy_face[iface,kk,ll] = ip_dest
                                        end
                                    end
                                end
                            end

                            for e=1:mesh.nelem
                                for ii=1:mesh.ngl
                                    for jj=1:mesh.ngl
                                        for kk=1:mesh.ngl
                                            ipp = mesh.connijk[e,ii,jj,kk]
                                            if (ipp > ip_kill)
                                                mesh.connijk[e,ii,jj,kk] -= 1
                                            end
                                            if (ipp == ip_kill)
                                                mesh.connijk[e,ii,jj,kk] = ip_dest
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
        ### work on single periodicity on the periodic2 boundaries if applicable
        if (size(per2_points,1) > 1)    
            for i=1:size(per2_points,1)
                found = false
                i1 = i+1
                while (i1 <= size(per2_points,1) && found == false)
                    ip = per2_points[i]
                    ip1 = per2_points[i1]
                    vec = [mesh.x[ip] - mesh.x[ip1], mesh.y[ip] - mesh.y[ip1], mesh.z[ip] - mesh.z[ip1]]
                    if (determine_colinearity(vec, nor2))
                        #@info "found a match", vec, nor2, ip, ip1,mesh.z[ip], mesh.z[ip1], mesh.y[ip], mesh.y[ip1], mesh.x[ip], mesh.x[ip1]
                        found = true
                        cond1 = mesh.x[ip]*abs(mesh.y[ip])*abs(mesh.z[ip]) < mesh.x[ip1]*abs(mesh.y[ip1])*abs(mesh.z[ip1])
                        cond2 = mesh.y[ip]*abs(mesh.x[ip])*abs(mesh.z[ip]) < mesh.y[ip1]*abs(mesh.x[ip1])*abs(mesh.z[ip1])
                        cond3 = mesh.z[ip]*abs(mesh.x[ip])*abs(mesh.y[ip]) < mesh.z[ip1]*abs(mesh.x[ip1])*abs(mesh.y[ip1])
                        if (cond1 || cond2 || cond3)
                            ip_dest = ip
                            ip_kill = ip1
                        else
                            ip_dest = ip1
                            ip_kill = ip
                        end
                        for e=1:mesh.nelem
                            for ii=1:mesh.ngl
                                for jj=1:mesh.ngl
                                    for kk=1:mesh.ngl
                                        ipp = mesh.connijk[e,ii,jj,kk]
                                        if (ipp == ip1 || ipp == ip)
                                            mesh.connijk[e,ii,jj,kk] = ip_dest
                                        end
                                    end
                                end
                            end
                        end
                        for iface = 1:size(mesh.poin_in_bdy_face,1)
                            for kk =1:mesh.ngl
                                for ll=1:mesh.ngl
                                    ipp = mesh.poin_in_bdy_face[iface,kk,ll]
                                    if (ipp == ip1 || ipp == ip)
                                        mesh.poin_in_bdy_face[iface,kk,ll] = ip_dest
                                    end
                                end
                            end
                        end
                        if !(ip_kill in mesh.connijk)
                            mesh.x[ip_dest] = min(mesh.x[ip_kill],mesh.x[ip_dest])
                            mesh.y[ip_dest] = min(mesh.y[ip_kill],mesh.y[ip_dest])
                            mesh.z[ip_dest] = min(mesh.z[ip_kill],mesh.z[ip_dest])
                            for ii=ip_kill:mesh.npoin-1
                                mesh.x[ii] = mesh.x[ii+1]
                                mesh.y[ii] = mesh.y[ii+1]
                                mesh.z[ii] = mesh.z[ii+1]
                            end
                            mesh.npoin = mesh.npoin-1
                            for iface =1:size(mesh.poin_in_bdy_face,1)
                                for kk=1:mesh.ngl
                                    for ll=1:mesh.ngl
                                        val = mesh.poin_in_bdy_face[iface,kk,ll]
                                        if (val > ip_kill)
                                            mesh.poin_in_bdy_face[iface,kk,ll] = val - 1
                                        end
                                        if (val == ip_kill)
                                            mesh.poin_in_bdy_face[iface,kk] = ip_dest
                                        end
                                    end
                                end
                            end

                            for e=1:mesh.nelem
                                for ii=1:mesh.ngl
                                    for jj=1:mesh.ngl
                                        for kk=1:mesh.ngl
                                            ipp = mesh.connijk[e,ii,jj,kk]
                                            if (ipp > ip_kill)
                                                mesh.connijk[e,ii,jj,kk] -= 1
                                            end
                                            if (ipp == ip_kill)
                                                mesh.connijk[e,ii,jj,kk] = ip_dest
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
        ### work on single periodicity on the periodic3 boundaries if applicable
        if (size(per3_points,1) > 1)
            for i=1:size(per3_points,1)
                found = false
                i1 = i+1
                while (i1 <= size(per3_points,1) && found == false)
                    ip = per3_points[i]
                    ip1 = per3_points[i1]
                    vec = [mesh.x[ip] - mesh.x[ip1], mesh.y[ip] - mesh.y[ip1], mesh.z[ip] - mesh.z[ip1]]
                    if (determine_colinearity(vec, nor3))
                        #@info "found a match", vec, nor3, ip, ip1, mesh.y[ip], mesh.y[ip1], mesh.x[ip], mesh.x[ip1], mesh.z[ip], mesh.z[ip1]
                        found = true
                        cond1 = mesh.x[ip]*abs(mesh.y[ip])*abs(mesh.z[ip]) < mesh.x[ip1]*abs(mesh.y[ip1])*abs(mesh.z[ip1])
                        cond2 = mesh.y[ip]*abs(mesh.x[ip])*abs(mesh.z[ip]) < mesh.y[ip1]*abs(mesh.x[ip1])*abs(mesh.z[ip1])
                        cond3 = mesh.z[ip]*abs(mesh.x[ip])*abs(mesh.y[ip]) < mesh.z[ip1]*abs(mesh.x[ip1])*abs(mesh.y[ip1])
                        if (cond1 || cond2 || cond3)
                            ip_dest = ip
                            ip_kill = ip1
                        else
                            ip_dest = ip1
                            ip_kill = ip
                        end
                        for e=1:mesh.nelem
                            for ii=1:mesh.ngl
                                for jj=1:mesh.ngl
                                    for kk=1:mesh.ngl
                                        ipp = mesh.connijk[e,ii,jj,kk]
                                        if (ipp == ip1 || ipp == ip)
                                            mesh.connijk[e,ii,jj,kk] = ip_dest
                                        end
                                    end
                                end
                            end
                        end
                        for iface = 1:size(mesh.poin_in_bdy_face,1)
                            for kk =1:mesh.ngl
                                for ll=1:mesh.ngl
                                    ipp = mesh.poin_in_bdy_face[iface,kk,ll]
                                    if (ipp == ip1 || ipp == ip)
                                        mesh.poin_in_bdy_face[iface,kk,ll] = ip_dest
                                    end
                                end
                            end
                        end
                        if !(ip_kill in mesh.connijk)
                            mesh.x[ip_dest] = min(mesh.x[ip_kill],mesh.x[ip_dest])
                            mesh.y[ip_dest] = min(mesh.y[ip_kill],mesh.y[ip_dest])
                            mesh.z[ip_dest] = min(mesh.z[ip_kill],mesh.z[ip_dest])
                            for ii=ip_kill:mesh.npoin-1
                                mesh.x[ii] = mesh.x[ii+1]
                                mesh.y[ii] = mesh.y[ii+1]
                                mesh.z[ii] = mesh.z[ii+1]
                            end
                            mesh.npoin = mesh.npoin-1
                            for iface =1:size(mesh.poin_in_bdy_face,1)
                                for kk=1:mesh.ngl
                                    for ll=1:mesh.ngl
                                        val = mesh.poin_in_bdy_face[iface,kk,ll]
                                        if (val > ip_kill)
                                            mesh.poin_in_bdy_face[iface,kk,ll] = val - 1
                                        end
                                        if (val == ip_kill)
                                            mesh.poin_in_bdy_face[iface,kk] = ip_dest
                                        end
                                    end
                                end
                            end

                            for e=1:mesh.nelem
                                for ii=1:mesh.ngl
                                    for jj=1:mesh.ngl
                                        for kk=1:mesh.ngl
                                            ipp = mesh.connijk[e,ii,jj,kk]
                                            if (ipp > ip_kill)
                                                mesh.connijk[e,ii,jj,kk] -= 1
                                            end
                                            if (ipp == ip_kill)
                                                mesh.connijk[e,ii,jj,kk] = ip_dest
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
        for e=1:mesh.nelem
            for i=1:mesh.ngl
                for j=1:mesh.ngl
                    for k=1:mesh.ngl
                        ip = mesh.connijk[e,i,j,k]
                    end
                end
            end
        end
    else
        #
        # 1D periodicity
        #
        if inputs[:AD] != FD()
            ip_dest = 1
            ip_kill = mesh.npoin_linear
            for e=1:mesh.nelem
                for i=1:mesh.ngl
                    if (mesh.connijk[e,i,1] == ip_kill)
                        mesh.connijk[e,i,1] = ip_dest
                    elseif (mesh.connijk[e,i,1] > ip_kill)
                        mesh.connijk[e,i,1] -= 1
                    end
                end
            end
            for ip=ip_kill:mesh.npoin-1
                mesh.x[ip] = mesh.x[ip+1]
            end
            mesh.npoin = mesh.npoin-1
        end
    end
    @info " periodicity_restructure!"
end
