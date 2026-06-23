function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)
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
            # Zero initial guess; exact field qe = manufactured solution u_ex
            # (used for the error/convergence check); 0 in the non-MMS modes.
            q.qn[ip,1] = 0.0
            q.qe[ip,1] = el_source_mode() == :mms ? manufactured_u(x, y) : 0.0
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
    # Zero initial guess; exact field = manufactured solution u_ex.
    qn[ip,1] = T(0.0)
    qe[ip,1] = T(MMS_A * sin(MMS_KX * xip) * cos(MMS_KY * yip))

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
                
                if x >= -0.75 && x <= 0.25 && y >= -0.75 && y <= 0.25 && old_ad_lvl[iel] < max_level
                    adapt_flags[iel] = refine_flag
                end
            end
        end
    end
    
end

function user_get_preadapt_flags!(adapt_flags, inputs, mesh, old_ad_lvl, connijk, nelem, ngl)

    max_level = inputs[:amr_max_level]
    for iel = 1:nelem
        m = 1
        for i = 1:ngl
            for j = 1:ngl
                ips = connijk[iel, i, j]
                
                # GEOMETRY HERE
                x = mesh.coords[ips, 1]
                y = mesh.coords[ips, 2]
                
                if x >= -0.9 && x <= -0.15 && y >= -0.95 && y <= -0.25 && old_ad_lvl[iel] < max_level
                    adapt_flags[iel] = refine_flag
                else
                    adapt_flags[iel] = nothing_flag
                end
            end
        end
    end
    
end
