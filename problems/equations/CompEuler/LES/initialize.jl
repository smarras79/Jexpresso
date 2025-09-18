using Random


function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    println_rank(" # Initialize fields for TCF with Ene equation ........................ "; msg_rank = rank)
    

    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    # 
    #---------------------------------------------------------------------------------
    qvars = ["ρ", "ρu", "ρv", "ρw", "ρe"]
    qoutvars = ["ρ", "u", "v", "w", "e"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------
    
    if (inputs[:backend] == CPU())
        PhysConst = PhysicalConst{Float64}()
        
        #
        # INITIAL STATE from scratch:
        #
        r0 = 1000.0 #m
        xc = 5000.0
        zc = 1500.0
        trc = 1.0
        new_param_set = create_updated_TD_Parameters(PhysConst.potential_temperature_reference_pressure)
        for ip = 1:mesh.npoin
            
            x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]
            r = sqrt( (x - xc)^2 + (z - zc)^2 )
            Δtr = 0.0 #K
            if r < r0
                Δtr = trc*(1.0 - r/r0)
            end
            
            u =10.0
            v = 0.0
            w = 0.0
            
            θ = 300.0
            p = 100000
            ρ = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p) #kg/m³

            θref = θ
            pref = p
            ρref = ρ
            
            q.qn[ip,1] = ρ
            q.qn[ip,2] = ρ*u
            q.qn[ip,3] = ρ*v
            q.qn[ip,4] = ρ*w
            q.qn[ip,5] = ρ*θ
            #q.qn[ip,6] = Δtr

            #Store initial background state for plotting and analysis of pertuebations
            q.qe[ip,1] = ρref
            q.qe[ip,2] = ρref*u
            q.qe[ip,3] = ρref*v
            q.qe[ip,4] = ρref*w
            q.qe[ip,5] = ρref*θref
            #q.qe[ip,6] = 0.0
        end
    end
    
    C = 0.3
    u = maximum(q.qn[:,2]./q.qn[:,1])
    v = maximum(q.qn[:,3]./q.qn[:,1])
    w = maximum(q.qn[:,4]./q.qn[:,1])
    speed = sqrt(u*u + v*v + w*w)
    inputs[:Δt] = min(C*mesh.Δelem_s/speed,  inputs[:Δt])
    println_rank(" Δt = ", inputs[:Δt]; msg_rank = rank)
    
    println_rank(" Initialize fields for TCF with Ene equation ........................ DONE "; msg_rank = rank)

    
    return q
end
