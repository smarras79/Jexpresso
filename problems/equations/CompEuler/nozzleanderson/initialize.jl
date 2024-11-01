function initialize(SD, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

    """
    @info " Initialize fields for 2D CompEuler with θ equation ........................ "
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    # 
    #---------------------------------------------------------------------------------
    qvars = ("ρ", "ρu", "E")       #These are the solution variables
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat; neqs=length(qvars))
    #---------------------------------------------------------------------------------

    γm1 = 0.4
    γ   = 1.4

    qout = copy(q.qn)
    ρ = 0.0
    T = 0.0

    lshock = false

    Ai = zeros(mesh.npoin)

    initial = zeros(mesh.npoin, 10)
    
    mass_flow = 0.59
    for iel_g = 1:mesh.nelem
        for i=1:mesh.ngl
            
            ip = mesh.connijk[iel_g,i,1]
            x  = mesh.x[ip]

            A = 1.0 + 2.2*(x - 1.5)^2
            initial[ip,1] = A
            initial[ip,5] = x
            if lshock
                if (x >= 0.0 && x < 0.5)
                    ρ = 1.0
                    T = 1.0
                elseif (x >= 0.5 && x < 1.5)
                    ρ = 1.0 - 0.366*(x - 0.5)
                    T = 1.0 - 0.167*(x - 0.5)
                elseif (x >= 1.5 && x < 2.1)
                    ρ = 0.634 - 0.702*(x - 1.5)
                    T = 0.833 - 0.4908*(x - 1.5)
                elseif (x >= 2.1 && x <= 3.0)
                    ρ = 0.5892 - 0.10228*(x - 1.5)
                    T = 0.93968 - 0.0622*(x - 1.5)
                end
            else
                if (x >= 0.0 && x <= 0.5)
                    ρ = 1.0
                    T = 1.0
                elseif (x > 0.5 && x <= 1.5)
                    ρ = 1.0 - 0.366*(x - 0.5)
                    T = 1.0 - 0.167*(x - 0.5)
                elseif (x > 1.5 && x <= 3.5)
                    ρ = 0.634 - 0.3879*(x - 1.5)
                    T = 0.833 - 0.3507*(x - 1.5)
                end
            end
            u = mass_flow/(ρ*A)
            e = T
            p = ρ*T
            
            initial[ip,2] = ρ
            initial[ip,3] = u
            initial[ip,4] = T
            
            
            q.qn[ip,1] = ρ*A
            q.qn[ip,2] = ρ*A*u
            q.qn[ip,3] = ρ*(e/γm1 + 0.5*γ*u*u)*A
            q.qn[ip,end] = p

            initial[ip,6] = q.qn[ip,1]
            initial[ip,7] = q.qn[ip,2]
            initial[ip,8] = q.qn[ip,3]
            

            #Store initial background state for plotting and analysis of pertuebations
            q.qe[ip,:] .= q.qn[ip,:]
            
            qout[ip,1] = q.qn[ip,1]/A
            qout[ip,2] = q.qn[ip,2]/q.qn[ip,1]
            qout[ip,3] = γm1*(q.qn[ip,3]/q.qn[ip,1] - 0.5*γ*u*u)
            qout[ip,end] = ρ*T
        end
    end

    initial[sortperm(initial[:, 5]), :]
    
    open("INITIAL.txt", "w") do io
        writedlm(io, round.([initial[:,5] initial[:,1] initial[:,2] initial[:,3] initial[:,4] initial[:,6] initial[:,7] initial[:,8]], digits=3))
    end
    
    varsout = ("ρ", "u", "e", "p")
    if (isempty(inputs[:outvars]) == false)
        varsout = inputs[:outvars]
    end
    for ivar=1:length(varsout)
        ivarname = varsout[ivar]
        plot_initial(SD, mesh.x, qout[:,ivar], ivarname, OUTPUT_DIR)
    end
    
    @info " Initialize fields for 2D CompEuler with θ equation ........................ DONE "
    
    return q
end
