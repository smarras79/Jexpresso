function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=3, ip=1)

    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1]
    ρu = q[2]
    ρE = q[3]

    u    = ρu/ρ
    E    = ρE/ρ
    Temp = (E - 0.5*u*u)/PhysConst.cv
    Press = perfectGasLaw_ρTtoP(PhysConst; ρ=ρ, Temp=Temp)

    F[1] = ρu
    F[2] = ρu*u + Press
    F[3] = (ρE + Press)*u
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)

    ρ  = q[1]
    ρu = q[2]
    ρE = q[3]

    u    = ρu/ρ
    E    = ρE/ρ
    Temp = (E - 0.5*u*u)/PhysConst.cv
    Press = ρ*PhysConst.Rair*Temp

    return T(ρu), T(ρu*u + Press), T((ρE + Press)*u)
end
