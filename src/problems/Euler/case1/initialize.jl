function initialize(SD::NSD_1D, PT::Euler, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

    """
    @info " Initialize fields for 1D Euler equations ........................ "
    
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, TFloat; neqs=3)

    #Cone properties:

    case = "sod"
    if (case === "sod")
        @info " Sod tube"

        ρL, PL, uL = 1.0,   1.0, 0.0
        ρR, PR, uR = 0.125, 0.1, 0.0
        
    	for iel_g = 1:mesh.nelem
            for i=1:mesh.ngl
                
                ip = mesh.conn[i,iel_g]
                x  = mesh.x[ip]

                if (x < 0.0)
                    ρ = ρL
                    P = PL
                    u = uL
                else
                    ρ = ρR
                    P = PR
                    u = uR
                end
                γ = 1.4
                q.qn[ip,1] = ρ                     #ρ
                q.qn[ip,2] = ρ*u                   #ρu
                q.qn[ip,3] = P/(γ - 1.0) + ρ*u*u/2 #ρE
                
            end
        end
    else
        error(" ERROR: Euler: initialize.jl: no initial conditions assigned"
    end
    

    @info "Initialize fields for system of 1D Euler equations ........................ DONE"

    return q
end

function initialize(SD::NSD_2D, PT::Euler, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
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
