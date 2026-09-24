#---------------------------------------------------------------------------------
# Inviscid fluxes of the 2D compressible Euler equations in total-energy
# form, q = (rho, rho u, rho v, rho E):
#
#   F = ( rho u,  rho u^2 + p,  rho u v,      (rho E + p) u )
#   G = ( rho v,  rho u v,      rho v^2 + p,  (rho E + p) v )
#
# with p = (gamma-1)(rho E - rho|u|^2/2).
#
# The viscous fluxes are NOT here: the deviatoric stress tensor and the
# heat flux are assembled by the kernel from the primitive variables that
# user_primitives.jl hands it (rhs.jl, _expansion_visc!(::NSD_2D)), with
# the coefficient coming from Sutherland's law plus DynSGS -- see the
# header of user_inputs.jl.
#
# Pressure is taken straight from the conservative state rather than
# through a temperature round-trip: across a shock the round-trip loses
# nothing analytically, but the direct form is one operation and cannot
# drift from the equation of state the DynSGS residual is built on.
#---------------------------------------------------------------------------------
function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρE = q[4]

    u = ρu/ρ
    v = ρv/ρ

    Pressure = PhysConst.γm1*(ρE - 0.5*(ρu*u + ρv*v))

    F[1] = ρu
    F[2] = ρu*u + Pressure
    F[3] = ρv*u
    F[4] = (ρE + Pressure)*u

    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v + Pressure
    G[4] = (ρE + Pressure)*v
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)

    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρE = q[4]

    u = ρu/ρ
    v = ρv/ρ

    Pressure = PhysConst.γm1*(ρE - T(0.5)*(ρu*u + ρv*v))

    return T(ρu), T(ρu*u + Pressure), T(ρv*u), T((ρE + Pressure)*u),
           T(ρv), T(ρu*v), T(ρv*v + Pressure), T((ρE + Pressure)*v)
end
