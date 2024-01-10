function periodicity_restructure!(mesh,inputs)

    #per1 = inputs[:per1dd]
    #per2 = inputs[:per2]
    #determine boundary vectors
    if (mesh.nsd == 2)  
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
            per1 = [0.0, 0.0]
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
            per2 = [0.0, 0.0]
        end     
        xx = zeros(size(mesh.x,1),1)
        yy = zeros(size(mesh.y,1),1)
        poin_bdy=zeros(Int64,size(mesh.poin_in_bdy_edge))
        xx .= mesh.x
        yy .= mesh.y
        poin_bdy .=mesh.poin_in_bdy_edge
        #=interval = [2,3,4]
        for iedge_bdy =1:size(mesh.bdy_edge_type,1)
        if (mesh.bdy_edge_type[iedge_bdy] == "periodic1" || mesh.bdy_edge_type[iedge_bdy] == "periodic2")
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
        for i=1:size(interval,1)
        if (interval[i] >= ip_kill)
        interval[i] -= 1
        end
        end
        end   
        end
        end
        end=# 
        # New periodicity interface
        for iedge_bdy =1:size(mesh.bdy_edge_type,1)
            if (mesh.bdy_edge_type[iedge_bdy] == "periodic1" || mesh.bdy_edge_type[iedge_bdy] == "periodic2")
                comp = mesh.bdy_edge_comp[iedge_bdy]
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
                            comp_per = mesh.bdy_edge_comp[iedge_per]
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
                                    if (ip_true < ip_true1)
                                        ip_dest = ip_true
                                        ip_kill = ip_true1
                                    elseif (ip_true >= ip_true1)
                                        ip_dest = ip_true1
                                        ip_kill = ip_true
                                    end
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
end
