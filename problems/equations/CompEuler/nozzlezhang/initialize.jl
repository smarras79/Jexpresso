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
    qvars = ("U1", "U2", "U3")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat; neqs=length(qvars))
    #---------------------------------------------------------------------------------

    PhysConst = PhysicalConst{Float64}()
    γ = 1.4
    cv = PhysConst.cv
    
    massflow = 0.579
    ρ0 = 1.2218
    p0 = 47892.4
    M0 = 1.5
    #T0 = 1.0 #these are the values
    u0 = M0*sqrt(γ*p0/ρ0);

    a = 0.5    
    for iel_g = 1:mesh.nelem
        for i=1:mesh.ngl
            
            ip = mesh.connijk[iel_g,i,1]
            x  = mesh.x[ip]
            A  = 1.0 + a*(x - 1.0)^2   # A/Athroat
            
            if (x >= -1.0 && x < -0.3)
                U1 = 1.0   #ρ
                U2 = 0.71  #ρu
                U3 = 2.752 #ρe
                
                q.qn[ip,1] = U1
                q.qn[ip,2] = U2
                q.qn[ip,3] = U3
            elseif (x >= -0.3 && x <= 1.0)
                U1 = 1.0 - 0.3*x        #ρ
                U2 = 0.71*(1.0 - 0.3*x) #ρu
                U3 = 2.752 - 0.2*x      #ρe
                
                q.qn[ip,1] = U1
                q.qn[ip,2] = U2
                q.qn[ip,3] = U3
            end
            q.qn[ip,end] = p0
            
            #Store initial background state for plotting and analysis of pertuebations
            q.qe[ip,:] .= q.qn[ip,:]
        end
    end
    
    for ivar=1:length(qvars)
        plot_initial(SD, mesh.x, q.qn[:,ivar], ivar, OUTPUT_DIR)
    end
    
    @info " Initialize fields for 1D adv diff ........................ DONE "
    
    return q
end
