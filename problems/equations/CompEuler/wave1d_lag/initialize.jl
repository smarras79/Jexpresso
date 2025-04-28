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
    qvars = ["u", "v"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))
    #---------------------------------------------------------------------------------
    if (inputs[:backend] == CPU()) 
        σ = Float64(0.15)
        σ2= σ*σ
        for iel_g = 1:mesh.nelem
            for i=1:mesh.ngl
            
                ip = mesh.connijk[iel_g,i,1]
                x = mesh.x[ip]

                ex = -(x)^2/σ2
                q.qn[ip,1] = 2^ex
                q.qn[ip,2] = 0.0

                #Store initial background state for plotting and analysis of pertuebations
                q.qe[ip,1] = 0.0
                q.qe[ip,2] = 0.0
            
            end
        end

        if (inputs[:lwrite_initial] == true)
            for ivar=1:length(qvars)
                plot_initial(SD, mesh.x, q.qn[:,ivar], ivar, OUTPUT_DIR)
            end
        end
    else
        σ = TFloat(0.15)
        σ2= σ*σ
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, mesh.x, σ2; ndrange = (mesh.npoin))
    end
    
    @info " Initialize fields for 1D adv diff ........................ DONE "
    
    return q
end

@kernel function initialize_gpu!(qn, qe, x, σ2)
    ip = @index(Global, Linear)

    T = eltype(x)
    xip = x[ip]
    
    ex = -(xip)^2/σ2
    qn[ip,1] = T(2^ex)
    qn[ip,2] = T(0.0)

    #Store initial background state for plotting and analysis of pertuebations
    qe[ip,1] = T(0.0)
    qe[ip,2] = T(0.0)

end
