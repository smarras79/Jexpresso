function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

            """
    @info " Initialize fields for 3D CompEuler with θ equation ........................ "
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    # 
    #---------------------------------------------------------------------------------
    qvars = ("u", "v")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))
    #---------------------------------------------------------------------------------
    
    if (inputs[:backend] == CPU())
        PhysConst = PhysicalConst{Float64}()
        if inputs[:lrestart] == true
            #
            # READ RESTART HDF5:
            #
            q.qn, q.qe = read_output(mesh.SD, inputs[:restart_input_file_path], inputs, mesh.npoin, HDF5(); nvar=length(qvars))
            PhysConst = PhysicalConst{Float64}()
        
            for ip=1:mesh.npoin
                ρ  = q.qn[ip,1]
                ρθ = q.qn[ip,end-1]
                θ  = ρθ/ρ
                P = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
                q.qn[ip,end] = P
            
                ρe  = q.qe[ip,1]
                ρθe = q.qe[ip,end-1]
                θe  = ρθe/ρ
                Pe = perfectGasLaw_ρθtoP(PhysConst, ρ=ρe, θ=θe)
                q.qe[ip,end] = Pe
            end
        
        else
            #
            # INITIAL STATE from scratch:
            #
            xc = 0.0#(maximum(mesh.x) + minimum(mesh.x))/2
            zc = 0.0
            yc = 0.0#(maximum(mesh.y) + minimum(me#m
            rc = 0.5 #m
        
            θc   =   1.0 #K
            for ip = 1:mesh.npoin
            
                x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]
            
                r = sqrt( (x - xc)^2 + (y - yc)^2 + (z - zc)^2 )
            
                Δθ = 0.0 #K
                if r < rc
                    Δθ = θc*(1.0 - r/rc)
                end
                q.qn[ip,1] = Δθ
                q.qn[ip,2] = 2*Δθ
                q.qe[ip,1] = 0.0
                q.qe[ip,2] = 0.0
                #end
            end
        end
    
    end
    return q
end

@kernel function initialize_gpu!(qn, qe, x, y, z, xc, rθ, zc, θref, θc, PhysConst, lpert)
    ip = @index(Global, Linear)

    T = eltype(x)
    x = x[ip]
    y = y[ip]
    z = z[ip]
    r = sqrt( (x - xc)^2 + (z - zc)^2 )
    Δθ = T(0.0) #K
    
    if r < rθ
        Δθ = T(θc*(T(1.0) - r/rθ))
    end
    
    θ = θref + Δθ
    p    = PhysConst.pref*(T(1.0) - PhysConst.g*z/(PhysConst.cp*θ))^(PhysConst.cpoverR) #Pa
    pref = PhysConst.pref*(T(1.0) - PhysConst.g*z/(PhysConst.cp*θref))^(PhysConst.cpoverR)
    ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p)    #kg/m³
    ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref) #kg/m³

    u = T(0.0)
    v = T(0.0)
    w = T(0.0)

    if (lpert)
        qn[ip,1] = ρ - ρref
        qn[ip,2] = ρ*u - ρref*u
        qn[ip,3] = ρ*v - ρref*v
        qn[ip,4] = ρ*w - ρref*w
        qn[ip,5] = ρ*θ - ρref*θref
        qn[ip,end] = p
    else
        qn[ip,1] = ρ
        qn[ip,2] = ρ*u
        qn[ip,3] = ρ*v
        qn[ip,4] = ρ*w
        qn[ip,5] = ρ*θ
        qn[ip,end] = p
    end

                    #Store initial background state for plotting and analysis of pertuebations
    qe[ip,1] = ρref
    qe[ip,2] = ρref*u
    qe[ip,3] = ρref*v
    qe[ip,4] = ρref*w
    qe[ip,5] = ρref*θref
    qe[ip,end] = pref

end
