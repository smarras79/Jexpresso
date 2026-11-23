function initialize(SD, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)
    """

    """
    println(" Initialize fields for 1D adv diff ........................ ")
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is used to allocate all necessary equation-dependent arrays
    # 
    #---------------------------------------------------------------------------------
    qvars = ("q")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))
    #---------------------------------------------------------------------------------
    
    if (inputs[:backend] == CPU())        
        for ip=1:mesh.npoin
             # centred so that q(x = 0) = exp(-26^2/18) ~ 5e-17: the inflow value
             # q = 0 then matches the initial state to machine precision
             z0 = 26
             sigma = 3.0
             q.qn[ip,1] =exp(-((mesh.coords[1,ip]-z0)^2)/(2*sigma^2))

            #Store initial background state for plotting and analysis of pertuebations
             q.qe[ip,1] = 0.0
            
        end

    else
        z0 = TFloat(26)
        sigma = TFloat(3.0)
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, @view(mesh.coords[1,:]), z0, sigma; ndrange = (mesh.npoin))
    end
    println(" Initialize fields for 1D adv diff ........................ DONE ")
    
    return q
end

@kernel function initialize_gpu!(qn, qe, x, z0, sigma)
    ip = @index(Global, Linear)

    T = eltype(x)
    xip = x[ip]

    qn[ip,1] =exp(-((xip-z0)^2)/(2*sigma^2))

            #Store initial background state for plotting and analysis of pertuebations
    qe[ip,1] = 0.0
end

