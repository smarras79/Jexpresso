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
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))
    #---------------------------------------------------------------------------------
    
    if (inputs[:backend] == CPU())        
        for ip=1:mesh.npoin
                
                q.qn[ip,1] = 0.0
                q.qn[ip,2] = 0.0

            #Store initial background state for plotting and analysis of pertuebations
                q.qe[ip,1] = 10.0
                q.qe[ip,2] = 0.0
            
        end

        if (inputs[:lwrite_initial] == true)
            for ivar=1:length(qvars)
                plot_initial(SD, mesh.x, q.qn[:,ivar], ivar, OUTPUT_DIR)
            end
        end
        
    else
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe; ndrange = (mesh.npoin))
    end
    @info " Initialize fields for 1D adv diff ........................ DONE "
    
    return q
end

@kernel function initialize_gpu!(qn, qe)
    ip = @index(Global, Linear)

    T=eltype(qn)
    qn[ip,1] = T(0.0)
    qn[ip,2] = T(0.0)
            #Store initial background state for plotting and analysis of pertuebations
    qe[ip,1] = T(10.0)
    qe[ip,2] = T(0.0)
end
