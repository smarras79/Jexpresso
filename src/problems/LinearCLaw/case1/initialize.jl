function initialize(SD::NSD_2D, PT::LinearCLaw, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """
    
    """
    @info " Initialize fields for system of Linear Conservation Laws ........................ "
    
    ngl   = mesh.ngl
    nsd   = mesh.nsd
    nelem = mesh.nelem
    npoin = mesh.npoin
    neqs  = 3

    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, neqs,TFloat)
    
    #Cone properties:
    c  = 1.0
    x0 = y0 = -0.8
    kx = ky = sqrt(2.0)/2.0
    ω  = 0.2
    d  = 0.5*ω/sqrt(log(2.0)); d2 = d*d
        
    for iel_g = 1:mesh.nelem
        for i=1:ngl
            for j=1:ngl

                ip = mesh.connijk[i,j,iel_g]
                x  = mesh.x[ip]
                y  = mesh.y[ip]
                p = q.qn[ip,1] = exp(- ((kx*(x - x0) + ky*(y - y0))^2)/d2) #p
                u = q.qn[ip,2] = kx*p/c                                    #u
                v = q.qn[ip,3] = ky*p/c                                    #v 
                
            end
        end
    end
    
    @info " Initialize fields for system of Linear Conservation Laws ........................ DONE"
    
    return q
end
