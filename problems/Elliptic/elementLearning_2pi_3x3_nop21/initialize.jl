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
            # Initial guess is zero. The exact field qe is the manufactured
            # solution (used for the error/convergence check against the
            # element-learning result); it is 0 in the non-MMS modes.
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
