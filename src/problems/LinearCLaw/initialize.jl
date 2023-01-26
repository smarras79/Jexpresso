include("../AbstractProblems.jl")

include("../../kernel/globalStructs.jl")
include("../../kernel/mesh/mesh.jl")
include("../../io/plotting/jeplots.jl")

function initialize(PT::LinearCLaw, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)

    @info " Initialize fields for system of Linear Conservation Laws ........................ "
    
    ngl   = mesh.nop + 1
    nsd   = mesh.nsd

    nvars = 3
    q = St_SolutionVars{TFloat}(zeros(mesh.npoin, 3),            # qn+1
                                zeros(mesh.npoin, 3),            # qn
                                zeros(1, 1),                     # qn-1
                                zeros(1, 1),                     # qn-2
                                zeros(1, 1),                     # qn-3
                                zeros(mesh.npoin, 3),            # qe
                                zeros(1, 1, 1, 1),               # qnelⁿ
                                zeros(1, 1, 1, 1),               # Fⁿ⁺¹
                                zeros(1, 1, 1, 1),               # Gⁿ⁺¹
                                zeros(1, 1, 1, 1))               # Hⁿ⁺¹
    
    #Cone properties:
    c  = 1.0
    x0 = y0 = -0.8
    kx = ky = sqrt(2)/2
    ω  = 0.2
    d  = 0.5*ω/sqrt(log(2)); d2 = d*d
        
    for iel_g = 1:mesh.nelem
        for i=1:ngl
            for j=1:ngl

                ip = mesh.connijk[i,j,iel_g]
                x  = mesh.x[ip]
                y  = mesh.y[ip]

                e = exp(- ((kx*(x - x0) + ky*(y - y0))^2)/d2)
                p = q.qn[ip,1] = e      #p
                u = q.qn[ip,2] = kx*e/c #u
                v = q.qn[ip,3] = ky*e/c #v
                
                # [ip] -> [i,j,iel]
                #q.F[i,j,iel_g,1] = c^2*u
                #q.F[i,j,iel_g,2] = p
                #q.F[i,j,iel_g,3] = 0
                #
                #q.G[i,j,iel_g,1] = c^2*v
                #q.G[i,j,iel_g,2] = 0
                #q.G[i,j,iel_g,3] = p
                
                #=e = exp(-log2((x*x + y*y)/0.06^2))
                p = q.qn[ip,1] = e      #p
                u = q.qn[ip,2] = 0.0 #u
                v = q.qn[ip,3] = 0.0 #v
                
                # [ip] -> [i,j,iel]
                q.F[i,j,iel_g,1] = c^2*u
                q.F[i,j,iel_g,2] = p
                q.F[i,j,iel_g,3] = 0

                q.G[i,j,iel_g,1] = c^2*v
                q.G[i,j,iel_g,2] = 0
                q.G[i,j,iel_g,3] = p
                =#
                
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
    for ivar=1:nvars
        title = string(varnames[ivar], ": initial condition")
        jcontour(mesh.x, mesh.y, q.qn[:,1], title, string(OUTPUT_DIR, "/", varnames[ivar], "-INIT.png"))
    end
    @info " Initialize fields for system of Linear Conservation Laws ........................ DONE"
    
    return q
end
