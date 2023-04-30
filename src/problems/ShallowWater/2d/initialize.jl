function initialize(SD::NSD_2D, PT::ShallowWater, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """
    
    """
    @info " Initialize fields for 2D Shallow water equations ........................ "
    
    ngl   = mesh.ngl
    nsd   = mesh.nsd
    nelem = mesh.nelem
    npoin = mesh.npoin
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, TFloat; neqs=3)
        
    #Cone properties:
       
    @info "Constant height and no flow shallow water" 
    for iel_g = 1:mesh.nelem
        for i=1:ngl
            for j=1:ngl

                ip = mesh.connijk[i,j,iel_g]
                x  = mesh.x[ip]
                y  = mesh.y[ip]
                q.qn[ip,1] = 0.75 + exp(-4*(abs(x)^2 + abs(y)^2))/4   #H
                q.qn[ip,2] = 0.001 * q.qn[ip,1]  #Hu
                q.qn[ip,3] = 0.001 * q.qn[ip,1]  #Hv 
                
            end
        end
    end
    
    @info "Initialize fields for system of Shallow Water equations ........................ DONE"
    
    return q
end

function initialize(SD::NSD_1D, PT::ShallowWater, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

    """
    @info "Initialize fields for system of 1D Shallow Water equations ........................ "

    ngl   = mesh.ngl
    nsd   = mesh.nsd
    nelem = mesh.nelem
    npoin = mesh.npoin
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, TFloat; neqs=3)

    #Cone properties:

    case = 1
    if (case == 1)
        @info "Constant height with a immersed bump SWASHES first steady state case"
    
    	for iel_g = 1:mesh.nelem
            for i=1:ngl
                for j=1:ngl

                    ip = mesh.conn[i,iel_g]
                    x  = mesh.x[ip]
                    Hb = bathymetry(x)
                    q.qn[ip,1] = 0.5                                    #H
                    q.qn[ip,2] = 0.0                                   #Hu
                end
            end
        end
    elseif (case == 2)
        @info "Constant height with a dryish bump SWASHES second steady state case"
    
        for iel_g = 1:mesh.nelem
            for i=1:ngl
                for j=1:ngl

                    ip = mesh.conn[i,iel_g]
                    x  = mesh.x[ip]
                    Hb = bathymetry(x)
                    q.qn[ip,1] = max(0.1,Hb)                                    #H
                    q.qn[ip,2] = 0.0                                   #Hu
                end
            end
        end
    elseif (case == 3)
        @info "SWASHES subcritical steady state case"
    
        for iel_g = 1:mesh.nelem
            for i=1:ngl
                for j=1:ngl

                    ip = mesh.conn[i,iel_g]
                    x  = mesh.x[ip]
                    Hb = bathymetry(x)
                    q.qn[ip,1] = 2.0
                    #if (x > mesh.xmin)                                    #H
                        q.qn[ip,2] = 0.0
                    #else
                     #   q.qn[ip,2] = 4.42
                    #end                                   #Hu
                end
            end
        end
    elseif (case == 4)
        @info "SWASHES transcritical no shock steady state case"

        for iel_g = 1:mesh.nelem
            for i=1:ngl
                for j=1:ngl

                    ip = mesh.conn[i,iel_g]
                    x  = mesh.x[ip]
                    Hb = bathymetry(x)
                    q.qn[ip,1] = 0.66
                    q.qn[ip,2] = 0.0
                end
            end
        end
    elseif (case == 5)
        @info "SWASHES transcritical with shock steady state case"

        for iel_g = 1:mesh.nelem
            for i=1:ngl
                for j=1:ngl

                    ip = mesh.conn[i,iel_g]
                    x  = mesh.x[ip]
                    Hb = bathymetry(x)
                    q.qn[ip,1] = 0.33
                    q.qn[ip,2] = 0.0
                end
            end
        end
    end
   
    @info "Initialize fields for system of 1D Shallow Water equations ........................ DONE"
    
    return q
end
