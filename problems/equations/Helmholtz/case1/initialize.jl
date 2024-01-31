function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

    """
    @info " Initialize fields for 2D Helmholtz equation ........................ "
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    # 
    #---------------------------------------------------------------------------------
    qvars = ("u")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat; neqs=length(qvars))
    #---------------------------------------------------------------------------------


        
    for ip =1:mesh.npoin
        x=mesh.x[ip]
        y=mesh.y[ip]           
        q.qn[ip,1] = sin(x/2)*exp(-x/2)*cos(y)

        q.qe[ip,1] = sin(x/2)*exp(-x/2)*cos(y)
    end
        
    
    outvarsref = ("u_ref")    

    @info " Initialize fields for 2D Helmholtz equation ........................ DONE "
    
    return q
end
