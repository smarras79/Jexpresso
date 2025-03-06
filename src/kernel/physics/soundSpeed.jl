function soundSpeed(npoin, integrator, SD)
    
    # Physical constants
    PhysConst = PhysicalConst{Float32}()
    pos::TInt = 2
    if (SD == NSD_2D())
        pos = 3
    elseif (SD == NSD_3D())
        pos = 4
    end
    # Initialize arrays
    ρ = integrator.u[1:npoin]
    θ = integrator.u[pos*npoin+1:(pos+1)*npoin]
    
    # Compute pressure using vectorized operation
    p = perfectGasLaw_ρθtoP(PhysConst, ρ, θ)
    
    # Compute speed of sound using vectorized operation
    c = sqrt.(PhysConst.γ .* p ./ ρ)
    
    # Find the maximum speed of sound
    max_c = maximum(c)
    
    return max_c
end

function computeCFL(npoin, neqs, dt, Δs, integrator, SD::NSD_1D; visc=(0.0))
    nothing
end

function computeCFL(npoin, neqs, dt, Δs, integrator, SD::NSD_2D; visc=(0.0))

    velmax = @SVector zeros(neqs)
    
    #u
    for ieq = 1:neqs-1
        idx = (ieq)*npoin
        velmax[i] = maximum(integrator.u[idx+1:ieq*npoin])
    end
    #u
#=    ieq = 2
    idx = (ieq-1)*npoin
    umax = maximum(integrator.u[idx+1:ieq*npoin])

    #umax = integrator.u[idx+1]
# Iterate over the specified range to find the maximum value
    #for i in (idx+2):(ieq*npoin)
    #    umax = max(umax, integrator.u[i])
    #end
    
    #v
    ieq = 3
    idx = (ieq-1)*npoin
    vmax = maximum(integrator.u[idx+1:ieq*npoin])=#  
    
    #velomax
    #velomax = max(umax, vmax)
    velomax = maximum(velmax)
    
    #speed of sound
    c     = soundSpeed(npoin, integrator, SD)
    
    cfl_u = velomax*dt/Δs #Advective CFL
    cfl_c = c*dt/Δs       #Acoustic CFL

    Δs2      = Δs*Δs
    μ        = maximum(visc)
    λ        = 2.0 #free parameter
    cfl_visc = dt*λ*μ/Δs2 #Viscous CFL

    println(" #  Advective CFL: ", cfl_u)
    println(" #  Acoustic  CFL: ", cfl_c)
    println(" #  Viscous   CFL: ", cfl_visc)
    
end

function computeCFL(npoin, neqs, dt, Δs, integrator, SD::NSD_3D; visc=(0.0))

    #u
    ieq = 2
    idx = (ieq-1)*npoin
    umax = maximum(integrator.u[idx+1:ieq*npoin])
    #v
    ieq = 3
    idx = (ieq-1)*npoin
    vmax = maximum(integrator.u[idx+1:ieq*npoin])
    #w
    ieq = 4
    idx = (ieq-1)*npoin
    wmax = maximum(integrator.u[idx+1:ieq*npoin])        
    
    #velomax
    velomax = max(umax, vmax, wmax)
    
    #speed of sound
    c     = soundSpeed(npoin, integrator, SD)
    
    cfl_u = velomax*dt/Δs #Advective CFL
    cfl_c = c*dt/Δs       #Acoustic CFL
    
    Δs2      = Δs*Δs
    μ        = maximum(visc)
    λ        = 2.0 #free parameter
    cfl_visc = dt*λ*μ/Δs2 #Viscous CFL
    
    println(" #  Advective CFL: ", cfl_u)
    println(" #  Acoustic  CFL: ", cfl_c)
    println(" #  Viscous   CFL: ", cfl_visc)
    
end
