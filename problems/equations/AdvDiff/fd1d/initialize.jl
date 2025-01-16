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
    qvars = ("q")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat; neqs=length(qvars))
    #---------------------------------------------------------------------------------
    
    σ = Float64(64.0)
    for iel_g = 1:mesh.nelem
        for i=1:mesh.ngl
            
            ip = mesh.connijk[iel_g,i,1]
            x = mesh.x[ip]
            
            q.qn[ip,1] = exp(-200.0*(x - 0.5)^2)

            #Store initial background state for plotting and analysis of pertuebations
            q.qe[ip,1] = 0.0
            
        end
    end

    varnames = ["q1"]
    write_output(NSD_1D(), q.qn, mesh, OUTPUT_DIR, inputs, varnames, PNG())
    
    @info " Initialize fields for 1D adv diff ........................ DONE "
    
    return q
end
