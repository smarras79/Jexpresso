include("../AbstractProblems.jl")

include("../../kernel/globalStructs.jl")
include("../../kernel/mesh/mesh.jl")
include("../../io/plotting/jeplots.jl")

function initialize(PT::LinearCLaw, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """
    
    """
    @info " Initialize fields for system of Linear Conservation Laws ........................ "
    
    ngl   = mesh.ngl
    nsd   = mesh.nsd
    nelem = mesh.nelem
    npoin = mesh.npoin
    neqs  = 3

    q = allocate_q(nelem, npoin, ngl, neqs)
    
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
                u = q.qn[ip,2] = kx*p/c                                      #u
                v = q.qn[ip,3] = ky*p/c                                      #v 
                
                
                #p = q.qn[ip,1] = exp(-(log(2))*((x*x + y*y)/0.06^2))
                #u = q.qn[ip,2] = 0.0
                #v = q.qn[ip,3] = 0.0
                
                # [ip] -> [i,j,iel]
                q.F[i,j,iel_g,1] = c^2*u
                q.F[i,j,iel_g,2] = p
                q.F[i,j,iel_g,3] = 0

                q.G[i,j,iel_g,1] = c^2*v
                q.G[i,j,iel_g,2] = 0
                q.G[i,j,iel_g,3] = p
                
            end
        end
    end
     
    #------------------------------------------
    # Plot initial condition:
    # Notice that I scatter the points to
    # avoid sorting the x and q which would be
    # becessary for a smooth curve plot.
    #------------------------------------------
    varnames = ["p", "u", "v"]
    for ivar=1:length(varnames)
        title = string(varnames[ivar], ": initial condition")
        jcontour(mesh.x, mesh.y, q.qn[:,ivar], title, string(OUTPUT_DIR, "/", varnames[ivar], "-INIT.png"))
    end
    @info " Initialize fields for system of Linear Conservation Laws ........................ DONE"
    
    return q
end
