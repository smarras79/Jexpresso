function initialize(SD, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

                        """
    @info " Initialize fields for 1D Sod tube ........................ "

    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is used to allocate all necessary equation-dependent arrays
    #
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρe"]
    qoutvars = ["ρ", "u", "p", "T"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------

    @info " Initialize fields for 1D Trixi wave  ........................ "

    PhysConst = PhysicalConst{Float64}()
    if (inputs[:backend] == CPU())


    	for ip = 1:mesh.npoin
        x  = mesh.x[ip]
        ρ, v, p = initial_wave(x)
        #ρ, v, p = initial_condition_weak_blast_wave(x)

	    q.qn[ip,1] = ρ
        q.qn[ip,2] = ρ*v
        q.qn[ip,3] = p/(PhysConst.γ - 1.0) + 0.5*ρ*v^2

        end
    else
        @mystop(" Error: no GPU yet for this case")
    end

    @info " Initialize fields for 1D Trixi wave  ........................ DONE "

    return q
end

# From Trixi
function initial_wave(x)
    if abs(x - 1) < 0.25
        ρ = 1.1691
        v = 0.1882 * sign(x - 1)
        p = 1.245
    else
        ρ = 1.0
        v = 0.0
        p = 1.0
    end
    return ρ, v, p
end

# From Trixi
function initial_condition_weak_blast_wave(x)
    # From Hennemann & Gassner JCP paper 2020 (Sec. 6.3)
    # Set up polar coordinates
    RealT = eltype(x)
    inicenter = SVector(0)
    x_norm = x[1] - inicenter[1]
    r = abs(x_norm)
    # The following code is equivalent to
    # phi = atan(0.0, x_norm)
    # cos_phi = cos(phi)
    # in 1D but faster
    cos_phi = x_norm > 0 ? 1 : -1

    # Calculate primitive variables
    rho = r > 0.5f0 ? one(RealT) : convert(RealT, 1.1691)
    v1 = r > 0.5f0 ? zero(RealT) : convert(RealT, 0.1882) * cos_phi
    p = r > 0.5f0 ? one(RealT) : convert(RealT, 1.245)

    return rho, v1, p
end
