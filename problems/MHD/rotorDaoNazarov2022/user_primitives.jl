#---------------------------------------------------------------------------------
# Conserved-form DynSGS (:dsgs_conserved): the kernel diffuses (ρ, ρu, ρv, E, ρw, B, ψ) themselves.
# With :dsgs_nazarov_energy the thermal part p/(γ−1) of E goes to slot neqs+2 (SGS.jl, dsgs_split_energy).
#---------------------------------------------------------------------------------
@inline rotor_energy_split(u) = (eth = pressure_mhd(u[1], u[2], u[3], u[5], u[4], u[6], u[7], u[8], u[9])/(γ_mhd - 1.0); (u[4] - eth, eth))

function user_primitives!(u, qe, uprimitive, ::TOTAL)
    for ieq = 1:9
        uprimitive[ieq] = u[ieq]
    end
    if dsgs_split_energy[]
        uprimitive[4], uprimitive[11] = rotor_energy_split(u)
    end
end

function user_primitives(u, qe, uprimitive, ::TOTAL)
    e4 = dsgs_split_energy[] ? rotor_energy_split(u)[1] : u[4]
    return SVector(u[1], u[2], u[3], e4, u[5], u[6], u[7], u[8], u[9])
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    error(" problems/MHD/rotorDaoNazarov2022: PERT() solution variables are not supported.")
end

# Output: ρ, u, v, w, p, Bx, By, Bz, ψ, T = p/ρ, pmag = ½|B|², Mach = |v|/√(γp/ρ) (Tóth 2000 plots ρ, p, pmag, Mach)
function user_uout!(ip, ET, uout, u, qe; kwargs...)
    ρ  = u[1]
    p  = pressure_mhd(ρ, u[2], u[3], u[5], u[4], u[6], u[7], u[8], u[9])
    vx, vy, vz = u[2]/ρ, u[3]/ρ, u[5]/ρ
    uout[1]  = ρ
    uout[2]  = vx
    uout[3]  = vy
    uout[4]  = vz
    uout[5]  = p
    uout[6]  = u[6]
    uout[7]  = u[7]
    uout[8]  = u[8]
    uout[9]  = u[9]
    uout[10] = p/ρ
    uout[11] = 0.5*(u[6]*u[6] + u[7]*u[7] + u[8]*u[8])
    uout[12] = sqrt(vx*vx + vy*vy + vz*vz)/sqrt(max(γ_mhd*p/ρ, 1e-300))
end
