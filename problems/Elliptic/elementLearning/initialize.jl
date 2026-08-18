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
            x=mesh.coords[1, ip]
            y=mesh.coords[2, ip]           
            q.qn[ip,1] = 0.0 #sin(x/2)*exp(-x/2)*cos(y)

            q.qe[ip,1] = 0.0 #sin(x/2)*exp(-x/2)*cos(y)
        end
    else
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, view(mesh.coords, 1, :), view(mesh.coords, 2, :); ndrange = mesh.npoin)
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
