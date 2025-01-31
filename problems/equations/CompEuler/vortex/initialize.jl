function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
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
    qvars = ("ρ", "ρu", "ρv", "ρe")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))
    #---------------------------------------------------------------------------------
    PhysConst = PhysicalConst{Float64}()

    
    mycase   = "vor"
    Rair     = PhysConst.Rair
    cv       = PhysConst.cv
    γ        = PhysConst.γ
    γinv     = 1.0/γ
    γm1      = γ - 1.0
    γm1inv   = 1.0/γm1
    γstar    = γ*γm1inv
    
    if inputs[:lrestart] == true
        #
        # READ RESTART HDF5:
        #
        q.qn, q.qe = read_output(mesh.SD, inputs[:restart_input_file_path], inputs, mesh.npoin, HDF5(); nvar=length(qvars))
        PhysConst = PhysicalConst{Float64}()
        for ip=1:mesh.npoin
            ρ  = q.qn[ip,1]
            ρθ = q.qn[ip,4]
            θ  = ρθ/ρ
            P = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
            q.qn[ip,end] = P
            
            ρe  = q.qe[ip,1]
            ρθe = q.qe[ip,4]
            θe  = ρθe/ρ
            Pe = perfectGasLaw_ρθtoP(PhysConst, ρ=ρe, θ=θe)
            q.qe[ip,end] = Pe
        end
        
    else
        if mycase=="vor"
            
            # base flow
            xc = 0.0
            yc = 0.0
            
            #
            # INITIAL STATE from scratch:
            #
            if inputs[:case] == "shu" #Shu
                
                αdeg = 0.0
                Minf = sqrt(2.0/γ)
                ρinf = 1.0
                pinf = 1.0
                Tinf = 1.0
                R    = 1.0
                σ    = 1.0
                β    = Minf*5.0*sqrt(2.0)/(4π)*exp(0.5)
                                
            elseif inputs[:case] == "vincent" #vincent
                
                αdeg = 90.0
                Minf = 0.4
                ρinf = 1.0
                pinf = 1.0/(γ*Minf*Minf)
                Tinf = 1.0
                R    = 1.5
                σ    = 1.0
                β    = Minf*27.0/(4π)*exp(2.0/9.0)
                
            elseif inputs[:case] == "hw" #Hesthaven and Warburton
                
                αdeg = 0.0
                Minf = sqrt(1.0/γ)
                ρinf = 1.0
                pinf = 1.0
                Tinf = 1.0
                R    = 0.7071067811865476
                σ    = 1.0
                β    = Minf*5.0/(2π)*exp(1.0)
                
            end
            
            α        = αdeg*π/180.0
            MinfCos  = Minf*cos(α)
            MinfSin  = Minf*sin(α)
            R2       = R*R
            σ2       = σ*σ            
            oneOverR = 1.0/R
            Const    = -0.5/(σ2*R2)
            g        = 1.0/(γ - 1.0)
            Tc       = 5.0
            
            for ip = 1:mesh.npoin
                
                x = mesh.x[ip]
                y = mesh.y[ip]
                
                f  = Const*( (x - xc)^2 +  (y - yc)^2 )
                Ω  = β*exp(f)
                
                du = -y*Ω*oneOverR
                dv =  x*Ω*oneOverR
                dT = -0.5*γm1*Ω^2

                ρ = (1.0 + dT)^γm1inv
                u = MinfCos + du
                v = MinfSin + dv                
                p = γinv*(1.0 + dT)^γstar
                
                T = p/(ρ*Rair)
                e = p/(γ - 1.0) + 0.5*ρ*(u^2 + v^2);
                
                q.qn[ip,1] = ρ
                q.qn[ip,2] = ρ*u
                q.qn[ip,3] = ρ*v
                q.qn[ip,4] = ρ*e
                q.qn[ip,end] = p
                
                #Store initial background state for plotting and analysis of pertuebations
                q.qe[ip,1] = 0.0
                q.qe[ip,2] = 0.0
                q.qe[ip,3] = 0.0
                q.qe[ip,4] = 0.0
                q.qe[ip,end] = 0.0
            end
            
        elseif mycase=="the"

            xc = (maximum(mesh.x) + minimum(mesh.x))/2
            yc = 0.0 #m
            r0 = 2.0 #m
            
            θref = 1.0 #K
            θc   = θref*0.01 #K
            ##
            for ip = 1:mesh.npoin
            
                x, y = mesh.x[ip], mesh.y[ip]
                r = sqrt( (x - xc)^2 + (y - yc)^2 )
            
                Δθ = 0.0 #K
                if r < r0
                    Δθ = θc*(1.0 - r/r0)
                end
                θ    = θref + Δθ
                p    = 1.0 #PhysConst.pref*(1.0 - PhysConst.g*y/(PhysConst.cp*θ))^(PhysConst.cpoverR) #Pa
                pref = 1.0 #PhysConst.pref*(1.0 - PhysConst.g*y/(PhysConst.cp*θref))^(PhysConst.cpoverR)
                ρ    = 1.0 #perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p)    #kg/m³
                ρref = 1.0 #perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref) #kg/m³
                
                u = 2.0
                v = 0.0

                q.qn[ip,1]   = ρ
                q.qn[ip,2]   = ρ*u
                q.qn[ip,3]   = ρ*v
                q.qn[ip,4]   = ρ*θ
                q.qn[ip,end] = p

                #Store initial background state for plotting and analysis of pertuebations
                q.qe[ip,1]   = ρref
                q.qe[ip,2]   = u
                q.qe[ip,3]   = v
                q.qe[ip,4]   = ρref*θref
                q.qe[ip,end] = pref
            end
            ##
        end
    end
    
    if (inputs[:lwrite_initial] == true)
        outvarsref = ("rho_init", "u_init", "v_init", "e_init", "p_init")
        write_vtk_ref(SD, mesh, q.qn.-q.qe, "initial_state", inputs[:output_dir]; nvar=length(q.qn[1,:]), outvarsref=outvarsref)
        
        outvarsref = ("rho_ref", "u_ref", "v_ref", "e_ref", "p_ref")    
        write_vtk_ref(SD, mesh, q.qe, "REFERENCE_state", inputs[:output_dir]; nvar=length(q.qe[1,:]), outvarsref=outvarsref)
    end
    
    @info " Initialize fields for 2D CompEuler with θ equation ........................ DONE "

    return q
end
