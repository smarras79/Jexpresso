function initialize(SD, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

                        """
    @info " Initialize fields for 1D Sod tube ........................ "
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is used to allocate all necessary equation-dependent arrays
    # 
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρe"]
    qoutvars = ["ρ", "u", "p", "T"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------
    
    @info " Initialize fields for 1D Trixi wave  ........................ "

    PhysConst = PhysicalConst{Float64}()
    if (inputs[:backend] == CPU()) 

        
    	for ip = 1:mesh.npoin
            
            x  = mesh.x[ip]

            if abs(x - 1) < 0.25
		ρ = 1.1691
		v = 0.1882 * sign(x - 1)
		p = 1.245
	    else
		ρ = 1.0
		v = 0.0
		p = 1.0
	    end
	    q.qn[ip,1] = ρ
            q.qn[ip,2] = ρ*v
            q.qn[ip,3] = p/(PhysConst.γ - 1.0) + 0.5*ρ*v^2
            
        end
    else
        @mystop(" Error: no GPU yet for this case")
    end

    @info " Initialize fields for 1D Trixi wave  ........................ DONE "

    return q
end
