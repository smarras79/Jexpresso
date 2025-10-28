function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

    """
    @info " Initialize fields for 2D incompressible NS of vorticity-stream function form ........................ "
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    # 
    #---------------------------------------------------------------------------------
    qvars = ("ω")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))
    #---------------------------------------------------------------------------------
    if (inputs[:backend] == CPU())    
        PhysConst = PhysicalConst{Float64}()
        if (inputs[:case] === "rtb")
        
            for iel_g = 1:mesh.nelem
                for j=1:mesh.ngl, i=1:mesh.ngl
                
                    ip = mesh.connijk[iel_g,i,j]
                    x, y = mesh.x[ip], mesh.y[ip]

                    if inputs[:SOL_VARS_TYPE] == PERT()
                        q.qn[ip] = cos(pi*x/2)*cos(pi*y/2)
                    
                        #Store initial background state for plotting and analysis of pertuebations
                        q.qe[ip] = 0

                    else
                        q.qn[ip] = cos(pi*x/2)*cos(pi*y/2)

                        #Store initial background state for plotting and analysis of pertuebations
                        q.qe[ip] = 0
                    end
                end
            end
        
        else
            error(" ERROR: CompEuler: initialize.jl:\n assign value to inputs[:case]")
        end

        #
        # Write reference to VTK:
        #  
        if (inputs[:lwrite_initial] == true)
            outvarsref = Array{Union{Nothing, String}}(nothing, q.neqs)
            for i = 1:length(outvarsref)
                outvarsref[i] = string(qvars[i], "_ref")
            end
            write_vtk_ref(SD, mesh, q.qe, "REFERENCE_state", inputs[:output_dir]; nvar=length(q.qe[1,:]), outvarsref=outvarsref)
        end
    
    end
    @info " Initialize fields for 2D incompressible NS of vorticity-stream function form ........................ DONE "
    
    return q
end

#=@kernel function initialize_gpu!(qn, qe, x, y, xc, rθ, yc, θref, θc, PhysConst,lpert)
    ip = @index(Global, Linear)

    T = eltype(x)
    xc2   = T(-2000.0)
    yc2   =  T(4000.0) #m
    r02   =  T(1000.0) #m
    qtrc2 =     T(1.0) #K

    θref = T(300.0) #K
    θc   =   T(2.0) #K
    qtrref = T(0.0) #K
    qtrc   = T(1.0) #K

    x = x[ip]
    y = y[ip]
    r = sqrt( (x - xc)^2 + (y - yc)^2 )
    Δθ = T(0.0) #K
    Δqtr = T(0.0)
    if r < rθ
        Δθ = T(θc*(T(1.0) - r/rθ))
        Δqtr = qtrc
    end

    r = sqrt( (x - xc2)^2 + (y - yc2)^2 )
    Δqtr2 = T(0.0) #K
    if r < r02
        Δqtr2 = qtrc2
    end
    qtr  = qtrref + Δqtr
    qtr2 = qtrref + Δqtr2

    θ = θref + Δθ
    p    = PhysConst.pref*(T(1.0) - PhysConst.g*y/(PhysConst.cp*θ))^(PhysConst.cpoverR) #Pa
    pref = PhysConst.pref*(T(1.0) - PhysConst.g*y/(PhysConst.cp*θref))^(PhysConst.cpoverR)
    ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p)    #kg/m³
    ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref) #kg/m³

    u = T(0.0)
    v = T(0.0)
    if (lpert)
        qn[ip,1] = ρ - ρref
        qn[ip,2] = ρ*u
        qn[ip,3] = ρ*v
        qn[ip,4] = ρ*θ - ρref*θref
        qn[ip,5] = qtr - qtrref
        qn[ip,6] = qtr2 - qtrref
        qn[ip,end] = p
    else
        qn[ip,1] = ρ
        qn[ip,2] = ρ*u
        qn[ip,3] = ρ*v
        qn[ip,4] = ρ*θ
        qn[ip,5] = qtr
        qn[ip,6] = qtr2 
        qn[ip,end] = p
    end

                    #Store initial background state for plotting and analysis of pertuebations
    qe[ip,1] = ρref
    qe[ip,2] = u
    qe[ip,3] = v
    qe[ip,4] = ρref*θref
    qe[ip,5] = qtrref
    qe[ip,6] = qtrref
    qe[ip,end] = pref

end=#
