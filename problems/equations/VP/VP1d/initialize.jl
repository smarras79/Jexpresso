function initialize(SD, PT, mesh::St_mesh, mesh_extra::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

    """
    @info " Initialize fields for 1D + 1D Vlasov-Poisson equations ........................ "
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is used to allocate all necessary equation-dependent arrays
    # 
    #---------------------------------------------------------------------------------
    qvars = ("u")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh_extra.npoin, mesh.ngl, qvars, TFloat, inputs[:backend], inputs[:l_extra_coord]; neqs=length(qvars))
    #---------------------------------------------------------------------------------
    if (inputs[:backend] == CPU()) 
       epsilon = 0.005
       beta = 0.2
       v_avg  = 2.4

        for iel = 1:mesh.nelem
            for i=1:mesh.ngl

                ip = mesh.connijk[iel,i,1]
                x = mesh.x[ip]

                for iel_extra = 1:mesh_extra.nelem
                    for j = 1:mesh_extra.ngl
            
                    ip_extra = mesh_extra.connijk[iel_extra,j,1]
                    v = mesh_extra.x[ip_extra]

                    ind = (ip-1)*mesh_extra.npoin + ip_extra

                    mu = ( exp(-(v-v_avg)^2/2) + exp(-(v+v_avg)^2/2) )/(2*sqrt(2*pi))
                    cs = cos(beta*x)

                    q.qn[ind,1] = epsilon*mu*cs

                    #Store initial background state for plotting and analysis of pertuebations
                    q.qe[ind,1] = 0.0

                    end            
                end
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
    
    @info " Initialize fields for 1D + 1D Vlasov-Poisson ........................ DONE "
    
    return q
end

