function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

        """
    @info " Initialize fields for 2D AD ........................ "
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    # 
    #---------------------------------------------------------------------------------
    qvars = ("q")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat; neqs=length(qvars))
    #---------------------------------------------------------------------------------
    xc = -0.25
    yc = 0.0
    s = 32.0
    for ip = 1:mesh.npoin
        
        x, y = mesh.x[ip], mesh.y[ip]
        
        a1 = ((x - xc))^2
        a2 = ((y - yc))^2
        q.qn[ip,1] = exp(-s*(a1 + a2))
        
        #Store initial background state for plotting and analysis of pertuebations
        q.qe[ip,1] = 0.0
    end
    
    #
    # Write reference and initial to VTK:
    #
    outvarsref = Array{Union{Nothing, String}}(nothing, q.neqs)
    for i = 1:length(outvarsref)
        outvarsref[i] = string(qvars[i], "_ref")
    end
    write_vtk_ref(SD, mesh, q.qe, "REFERENCE_state", inputs[:output_dir]; nvar=length(q.qe[1,:]), outvarsref=outvarsref)

    for i = 1:length(outvarsref)
        outvarsref[i] = string(qvars[i], "_init")
    end
    write_vtk_ref(SD, mesh, q.qn, "initial_state", inputs[:output_dir]; nvar=length(q.qn[1,:]), outvarsref=outvarsref)
    
    
    @info " Initialize fields for 2D AD ........................ DONE "

    return q
end
