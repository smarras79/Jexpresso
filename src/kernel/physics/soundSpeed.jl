function soundSpeed(npoin, integrator)
    
    #speed of sound
    PhysConst = PhysicalConst{TFloat}()
    ctmp = Float32(0.0)
    c = Float32(0.0)
    for i=1:npoin
        ρ  = integrator.u[i]
        θ  = integrator.u[3*npoin+i]
        p  = perfectGasLaw_ρθtoP(PhysConst; ρ=ρ, θ=θ)
        c  = sqrt(PhysConst.γ*p/ρ)
        c  = max(c, ctmp)
        ctmp = c
    end

    return c
    
end

function computeCFL(npoin, dt, Δs, integrator, SD::NSD_2D)

    #u
    ieq = 2
    idx = (ieq-1)*npoin
    umax = maximum(integrator.u[idx+1:ieq*npoin])
    #v
    ieq = 3
    idx = (ieq-1)*npoin
    vmax = maximum(integrator.u[idx+1:ieq*npoin])        
    
    #velomax
    velomax = max(umax, vmax)
    
    #speed of sound
    c     = soundSpeed(npoin, integrator)
    
    cfl_u = velomax*dt/Δs #Advective CFL
    cfl_c = c*dt/Δs       #Acoustic CFL

    println(" #  Advective CFL: ", cfl_u)
    println(" #  Acoustic  CFL: ", cfl_c)
    
end

function computeCFL(npoin, dt, Δs, integrator, SD::NSD_3D)

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
    vmax = maximum(integrator.u[idx+1:ieq*npoin])        
    
    #velomax
    velomax = max(umax, vmax, wmax)
    
    #speed of sound
    c     = soundSpeed(npoin, integrator)
    
    cfl_u = velomax*dt/Δs #Advective CFL
    cfl_c = c*dt/Δs       #Acoustic CFL

    println(" #  Advective CFL: ", cfl_u)
    println(" #  Acoustic  CFL: ", cfl_c)
    
end
