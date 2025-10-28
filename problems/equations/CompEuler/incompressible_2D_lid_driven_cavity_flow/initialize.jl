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
                        q.qn[ip] = 0
                        #Store initial background state for plotting and analysis of pertuebations
                        q.qe[ip] = 0

                    else
                        q.qn[ip] = 0
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



