function initialize(SD, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

    """
    @info " Initialize fields for 1D adv diff ........................ "
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is used to allocate all necessary equation-dependent arrays
    # 
    #---------------------------------------------------------------------------------
    qvars = ("h","u")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat; neqs=length(qvars))
    #---------------------------------------------------------------------------------
    
            
    for ip=1:mesh.npoin
            
            q.qn[ip,1] = 0.0
            q.qn[ip,2] = 0.0

            #Store initial background state for plotting and analysis of pertuebations
            q.qe[ip,1] = 10.0
            q.qe[ip,2] = 0.0
            
    end

    varnames = ["h","u"]
    write_output(NSD_1D(), q.qn, mesh, OUTPUT_DIR, inputs, varnames, PNG())
    
    @info " Initialize fields for 1D adv diff ........................ DONE "
    
    return q
end
