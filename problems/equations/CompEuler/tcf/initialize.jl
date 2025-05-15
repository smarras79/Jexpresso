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
        new_param_set = create_updated_TD_Parameters(Float64(101325.0))
        for ip = 1:mesh.npoin
            
            x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]
                      
            u = 2.0
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
            q.qn[ip,end] = p

            #Store initial background state for plotting and analysis of pertuebations
            q.qe[ip,1] = ρref
            q.qe[ip,2] = ρref*u
            q.qe[ip,3] = ρref*v
            q.qe[ip,4] = ρref*w
            q.qe[ip,5] = ρref*θref
            q.qe[ip,end] = pref
        end
    end
    
    println_rank(" Initialize fields for TCF with Ene equation ........................ DONE "; msg_rank = rank)

    return q
end
