function user_source!(S,
                      q,
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=1, x=0.0, y=0.0, ymin=0.0, ymax=1000.0, xmin=0.0, xmax=1000.0)

    PhysConst = PhysicalConst{Float64}()

    #
    # S(q(x)) = -ρg
    #
    ρ = q[1]

    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0

    #--------------
    # SPONGE (weak, top of domain, thickness = ymax - zsponge)
    # Adapted to 2D from problems/CompEuler/LESICP2/user_source.jl
    # where the sponge acts along the vertical (z in 3D, y here).
    #--------------
    if inputs[:lsponge] == true
        ys = inputs[:zsponge]
        α  = 0.05   # weak sponge strength
        if y >= ys
            betay_coe = α*sinpi(0.5*(y - ys)/(ymax - ys))
        else
            betay_coe = 0.0
        end
        cs = betay_coe

        S[2] -= cs*(q[2] - qe[2])
        S[3] -= cs*(q[3] - qe[3])
        S[4] -= cs*(q[4] - qe[4])
    end

end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=1, x=0.0, y=0.0, ymin=0.0, ymax=1000.0, xmin=0.0, xmax=1000.0)

    PhysConst = PhysicalConst{Float64}()

    #
    # S(q(x)) = -ρg
    #
    ρ = q[1]

    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0

    if inputs[:lsponge] == true
        ys = inputs[:zsponge]
        α  = 0.05
        if y >= ys
            betay_coe = α*sinpi(0.5*(y - ys)/(ymax - ys))
        else
            betay_coe = 0.0
        end
        cs = betay_coe

        S[2] -= cs*(q[2] - qe[2])
        S[3] -= cs*(q[3] - qe[3])
        S[4] -= cs*(q[4] - qe[4])
    end

end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::NCL,
                      ::AbstractPert;
                      neqs=1, x=0.0, y=0.0, ymin=0.0, ymax=1000.0, xmin=0.0, xmax=1000.0)

    PhysConst = PhysicalConst{Float64}()

    S[1] = 0.0
    S[2] = 0.0
    S[3] = -PhysConst.g
    S[4] = 0.0

end

function user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)

    T = eltype(q)
    ρ = q[1]
    return T(0.0), T(0.0), T(-ρ*PhysConst.g), T(0.0)
end
