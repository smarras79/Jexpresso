function initialize(SD, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)
    """
       Initial condition for the 1D inviscid Burgers equation:

           q(x,0) = sin(2π x / L),    L = xmax - xmin

        which steepens into a stationary shock at x = 0.5 at t = 1/(2π).
    """
    println(" Initialize fields for 1D Burgers + DSGS ................... ")

    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is used to allocate all necessary equation-dependent arrays
    #
    #---------------------------------------------------------------------------------
    qvars = ["q"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))
    #---------------------------------------------------------------------------------

    xmin = inputs[:xmin]
    xmax = inputs[:xmax]
    L    = xmax - xmin

    for iel_g = 1:mesh.nelem
        for i = 1:mesh.ngl

            ip = mesh.connijk[iel_g,i,1]
            x  = mesh.coords[1,ip]

            q.qn[ip,1] = sin(2.0*π*(x - xmin)/L)

            #Store initial background state for plotting and analysis of perturbations
            q.qe[ip,1] = 0.0

        end
    end

    println(" Initialize fields for 1D Burgers + DSGS ................... DONE ")

    return q
end
