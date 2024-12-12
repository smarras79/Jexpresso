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
    qvars = ("ρ", "u", "e")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))
    #---------------------------------------------------------------------------------

    case = "sod" # inputs[:case]
    
    PhysConst = PhysicalConst{Float64}()
    if (inputs[:backend] == CPU()) 

        if case == "sod"
            @info " Sod tube"

            ρL, uL, pL = 1.000, 0.0, 1.0
            ρR, uR, pR = 0.125, 0.0, 0.1
            xshock_initial = 0.5
            
    	    for iel_g = 1:mesh.nelem
                for i=1:mesh.ngl
                    
                    ip = mesh.connijk[iel_g,i,1]
                    x = mesh.x[ip]
                    
                    if (x < xshock_initial)
                        ρ = ρL
                        u = uL
                        p = pL
                    else
                        ρ = ρR
                        u = uR
                        p = pR
                    end
                    q.qn[ip,1] = ρ                       #ρ
                    q.qn[ip,2] = ρ*u                     #ρu
                    q.qn[ip,3] = p/(PhysConst.γ - 1.0) + 0.5*ρ*u*u #ρE
                    
                end
            end
            
        elseif case == "sound"
            xs  = 1.5
            u   = 0.0
            ωsq = 0.125^2
            for iel_g = 1:mesh.nelem
                for i=1:mesh.ngl
                    
                    ip = mesh.connijk[iel_g,i,1]
                    x = mesh.x[ip]

                    ρ = 1.0
                    p = exp(-log(2) * ((x - xs)^2)/ωsq) + 1.0
                    u = 0.0
                    
                    q.qn[ip,1] = ρ
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,2] = p/(PhysConst.γ - 1.0) + 0.5*ρ*u*u #ρE

                    q.qe[ip,:] = 0.0
                    
                end
            end
            
            if (inputs[:lwrite_initial] == true)
                for ivar=1:length(qvars)
                    plot_initial(SD, mesh.x, q.qn[:,ivar], ivar, OUTPUT_DIR)
                end
            end

        elseif (case === "sod")
            @info " Sod tube"

            ρL, uL, pL = 1.000, 0.0, 1.0
            ρR, uR, pR = 0.125, 0.0, 0.1
            xshock_initial = 0.5
            
    	    for iel_g = 1:mesh.nelem
                for i=1:mesh.ngl
                    
                    ip = mesh.connijk[i,iel_g]
                    x  = mesh.x[ip]
                    
                    if (x < xshock_initial)
                        ρ = ρL
                        u = uL
                        p = pL
                    else
                        ρ = ρR
                        u = uR
                        p = pR
                    end
                    q.qn[ip,1] = ρ                       #ρ
                    q.qn[ip,2] = ρ*u                     #ρu
                    q.qn[ip,3] = p/(PhysConst.γ - 1.0) + 0.5*ρ*u*u #ρE
                    
                end
            end
            
        end #case
    else
        @error(" SOD: I did not set up the GPU initial, bc, etc. for this case.")
        nothing
    end
    
    @info " Initialize fields for 1D adv diff ........................ DONE "
    
    return q
end

@kernel function initialize_gpu!(qn, qe, x, σ2)
    ip = @index(Global, Linear)

    T = eltype(x)
    xip = x[ip]
    
    ex = -(xip - 1)^2/σ2
    qn[ip,1] = T(2^ex)
    qn[ip,2] = T(0.0)

    #Store initial background state for plotting and analysis of pertuebations
    qe[ip,1] = T(0.0)
    qe[ip,2] = T(0.0)

end
