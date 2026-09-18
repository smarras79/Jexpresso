#---------------------------------------------------------------------------------
# Primitive variables handed to the viscous assembly, and the output map.
#
# uprimitive is what _expansion_visc!(::NSD_2D) differentiates.  On the
# total-energy path that kernel is the REAL compressible Navier-Stokes
# viscous operator: it builds the deviatoric stress
#
#     tau_ij = mu ( du_i/dx_j + du_j/dx_i - (2/3) delta_ij div u )
#
# from slots 2 and 3, adds the viscous work tau.u to the energy slot, and
# puts the conduction flux on the gradient of slot 4.  So the primitive set
# must be (rho, u, v, T) with T the temperature variable of Nazarov &
# Hoffman,
#
#     T = E/rho - |u|^2/2 = p/((gamma-1) rho),
#
# i.e. the SPECIFIC INTERNAL ENERGY -- the paper scales cv = 1, and in SI
# that quantity is cv*T_physical.  Putting the physical temperature in K
# here instead would silently rescale every slot-4 coefficient by
# cv = 717.5 and over-diffuse the energy equation by that factor.
#
# The two coefficients that ride on this convention:
#   DynSGS       kappa = Pr_a/(gamma-1) * mu      (:Pr, the artificial Pr)
#   Sutherland   k/cv  = gamma mu(T)/Pr           (:Pr_lam, the physical Pr)
# Both are applied by the kernel; see the header of user_inputs.jl.
#
# user_uout! is independent of all this: it writes the human-readable
# fields named by qoutvars in initialize.jl, where T IS the temperature in
# kelvin.
#---------------------------------------------------------------------------------
function user_primitives!(u, qe, uprimitive, ::TOTAL)

    ρ  = u[1]
    ρu = u[2]
    ρv = u[3]
    ρE = u[4]

    uprimitive[1] = ρ
    uprimitive[2] = ρu/ρ
    uprimitive[3] = ρv/ρ
    uprimitive[4] = ρE/ρ - 0.5*(ρu*ρu + ρv*ρv)/(ρ*ρ)   # T = p/((gamma-1) rho) = cv*T[K]
end

function user_primitives!(u, qe, uprimitive, ::PERT)

    ρ  = u[1] + qe[1]
    ρu = u[2] + qe[2]
    ρv = u[3] + qe[3]
    ρE = u[4] + qe[4]

    uprimitive[1] = ρ
    uprimitive[2] = ρu/ρ
    uprimitive[3] = ρv/ρ
    uprimitive[4] = ρE/ρ - 0.5*(ρu*ρu + ρv*ρv)/(ρ*ρ)
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)

    ρ  = u[1]
    ρu = u[2]
    ρv = u[3]
    ρE = u[4]

    return T(ρ), T(ρu/ρ), T(ρv/ρ), T(ρE/ρ - T(0.5)*(ρu*ρu + ρv*ρv)/(ρ*ρ))
end

#
# Output fields, in the order of qoutvars = ["rho", "u", "v", "p", "T"].
#
# T in kelvin is the field to look at on the wall: the paper's figure 2(a)
# and figure 3(c) plot the Stanton number, which is the wall-normal
# gradient of exactly this, and the cold-wall ratio T_w/T_inf = 2.34 is
# what sets the separation-bubble length.
#
function user_uout!(ip, ET, uout, u, qe; kwargs...)

    PhysConst = PhysicalConst{Float64}()

    ρ  = u[1]
    ρu = u[2]
    ρv = u[3]
    ρE = u[4]

    uout[1] = ρ                                              # rho
    uout[2] = ρu/ρ                                           # u
    uout[3] = ρv/ρ                                           # v
    uout[4] = PhysConst.γm1*(ρE - 0.5*(ρu*ρu + ρv*ρv)/ρ)     # p
    uout[5] = uout[4]/(PhysConst.Rair*ρ)                     # T [K]
end
