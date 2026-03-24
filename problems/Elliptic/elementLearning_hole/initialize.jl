function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

    """
    println(" # Initialize fields for 2D Helmholtz equation ........................ ")
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    # 
    #---------------------------------------------------------------------------------
    qvars = ["u"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))
    #---------------------------------------------------------------------------------


    if (inputs[:backend] == CPU())        
        for ip =1:mesh.npoin
            x=mesh.x[ip]
            y=mesh.y[ip]           
            q.qn[ip,1] = 0.0 #sin(x/2)*exp(-x/2)*cos(y)

            q.qe[ip,1] = 0.0 #sin(x/2)*exp(-x/2)*cos(y)
        end
    else
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, mesh.x, mesh.y; ndrange = mesh.npoin)
    end
        
    
    outvarsref = ("u_ref")    

     println(" # Initialize fields for 2D Helmholtz equation ........................ DONE ")
    
    return q
end

@kernel function initialize_gpu!(qn, qe, x, y)
    T = eltype(qn)

    ip = @index(Global, Linear)
    xip = x[ip]
    yip = y[ip]
    qn[ip,1] = sin(xip/2)*exp(-xip/2)*cos(yip)

    qe[ip,1] = sin(xip/2)*exp(-xip/2)*cos(yip)
        
end

function user_get_adapt_flags!(adapt_flags, inputs, mesh, old_ad_lvl, q, qe, connijk, nelem, ngl)

    max_level = inputs[:amr_max_level]
    for iel = 1:nelem
        m = 1
        for i = 1:ngl
            for j = 1:ngl
                ips = connijk[iel, i, j]
                
                # GEOMETRY HERE
                x = mesh.coords[ips, 1]
                y = mesh.coords[ips, 2]
                
                if x >= -0.75 && x =< 0.25 && y >= -0.75 && y =< 0.25 && old_ad_lvl[iel] < max_level
                    adapt_flags[iel] = refine_flag
                end
            end
        end
    end
    
end
