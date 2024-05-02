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
    qvars = ("u", "v")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat; neqs=length(qvars))
    #---------------------------------------------------------------------------------
    
    σ = Float64(0.15)
    σ2= σ*σ
    for iel_g = 1:mesh.nelem
        for i=1:mesh.ngl
            
            ip = mesh.connijk[iel_g,i,1]
            x = mesh.x[ip]

            ex = -(x - 1)^2/σ2
            q.qn[ip,1] = 2^ex
            q.qn[ip,2] = 0.0

            #Store initial background state for plotting and analysis of pertuebations
            q.qe[ip,1] = 0.0
            q.qe[ip,2] = 0.0
            
        end
    end
    
    for ivar=1:length(qvars)
        plot_initial(SD, mesh.x, q.qn[:,ivar], ivar, OUTPUT_DIR)
    end
    
    @info " Initialize fields for 1D adv diff ........................ DONE "
    
    return q
end
