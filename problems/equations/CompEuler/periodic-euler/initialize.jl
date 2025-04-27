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
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))
    #---------------------------------------------------------------------------------
    if (inputs[:backend] == CPU())
        xc = (mesh.xmax + mesh.xmin)/2
        yc = 8.0 #m
        sx = 1.0
        sy = 1.0
        A = 1.0
        for iel_g = 1:mesh.nelem
            for j=1:mesh.ngl, i=1:mesh.ngl
            
                ip = mesh.connijk[iel_g,i,j]
                x, y = mesh.x[ip], mesh.y[ip]
            
                a1 = -((x - xc)/sx)^2
                a2 = -((y - yc)/sy)^2
                q.qn[ip,1] = A*exp(a1)*exp(a2)
            
             #Store initial background state for plotting and analysis of pertuebations
                q.qe[ip,1] = 0.0
            end
        end
    

        #
        # Write reference to VTK:
        #
        outvarsref = Array{Union{Nothing, String}}(nothing, q.neqs)
        for i = 1:length(outvarsref)
            outvarsref[i] = string(qvars[i], "_ref")
        end
        write_vtk_ref(SD, mesh, q.qe, "REFERENCE_state", inputs[:output_dir]; nvar=length(q.qe[1,:]), outvarsref=outvarsref)
    else
        xc = TFloat((maximum(mesh.x) + minimum(mesh.x))/2)
        yc = TFloat(8.0)
        sx = TFloat(1.0)
        sy = TFloat(1.0)
        A = TFloat(1.0)
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, mesh.x, mesh.y, xc, yc, sx, sy, A; ndrange = (mesh.npoin))
    end
    @info " Initialize fields for 2D AD ........................ DONE "

    return q
end

@kernel function initialize_gpu!(qn, qe, x, y, xc, yc, sx, sy, A)
    ip = @index(Global, Linear)

    T = eltype(x)
    x = x[ip]
    y = y[ip]
    

    a1 = -((x - xc)/sx)^2
    a2 = -((y - yc)/sy)^2
    qn[ip,1] = A*exp(a1)*exp(a2)

             #Store initial background state for plotting and analysis of pertuebations
    qe[ip,1] = 0.0

end
