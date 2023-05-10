function initialize(SD::NSD_1D, ET::AdvDiff, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)

    println(" # Initialize fields for AdvDiff ........................")
    
    qinit = Array{TFloat}(undef, mesh.npoin, 1)
    q     = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, TFloat; neqs=1)
    
    σ = Float64(64.0)
    for iel_g = 1:mesh.nelem
        for i=1:mesh.ngl
            
            ip = mesh.connijk[i,iel_g]
            x  = mesh.x[ip]
            
            #q.qn[ip, 1] = exp(-σ*x*x)
            q.qn[ip, 1] = exp(-200.0*(x - 0.25)^2)
            
        end
    end
    
    #------------------------------------------
    # Plot initial condition:
    # Notice that I scatter the points to
    # avoid sorting the x and q which would be
    # becessary for a smooth curve plot.
    #------------------------------------------
    title = string( "Tracer: initial condition")
    plot_curve(mesh.x, q.qn[:,1], title, string(OUTPUT_DIR, "/INIT.png"))
    
    println(" # Initialize fields for AdvDiff ........................ DONE")
    
    return q
end


function initialize(SD::NSD_2D, ET::AdvDiff, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)

    println(" # Initialize fields for AdvDiff ........................")
        
    ngl  = mesh.nop + 1
    nsd  = mesh.nsd   
    q    = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, TFloat; neqs=1)
    
    test_case = "kopriva.5.3.5"
    #test_case = "giraldo.15.8"
    if (test_case == "kopriva.5.3.5")
        #Cone properties:
        ν = inputs[:νx] 
        if ν == 0.0
            ν = 0.01
        end
        σ = 1.0/ν
        (xc, yc) = (-0.5, -0.5)
        
        for iel_g = 1:mesh.nelem
            for i=1:ngl
                for j=1:ngl

                    ip = mesh.connijk[i,j,iel_g]
                    x  = mesh.x[ip]
                    y  = mesh.y[ip]

                    q.qn[ip,1] = exp(-σ*((x - xc)*(x - xc) + (y - yc)*(y - yc)))
                    u          = 0.8 #constant
                    v          = 0.8 #constant

                end
            end
        end
    elseif (test_case == "giraldo.15.8")
        σ = 32.0
        (xc, yc) = (-0.5, 0.0)
                
        for iel_g = 1:mesh.nelem
            for i=1:ngl
                for j=1:ngl

                    ip = mesh.connijk[i,j,iel_g]
                    x, y = mesh.x[ip], mesh.y[ip]

                    q.qn[ip,1] = exp(-σ*((x - xc)*(x - xc) + (y - yc)*(y - yc)))
                end
            end
        end
    end
    println(" # Initialize fields for AdvDiff ........................ DONE")
    
    return q
end
