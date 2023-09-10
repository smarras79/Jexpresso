using LinearSolve

# ASCII
function solution_norms(sol::ODESolution, OUTPUT_DIR::String, inputs::Dict;)
    
    println(string(" # Writing diagnostics to ASCII file:", OUTPUT_DIR, "*.dat ...  ") )

    fname = @sprintf "diagnostics.dat"
    open(string(OUTPUT_DIR, "/", fname), "w") do f
        @printf(f, "time |  norm₁ | norm₂ | norm∞ \n")
        for iout = 1:inputs[:ndiagnostics_outputs]
            #@printf(f, " %.6f %.6f %.6f %.6f \n", sol.t[iout], norm(sol.u[iout][:], 1), norm(sol.u[iout][:], 2), norm(sol.u[iout][:], Inf))
            @printf(f, " %d %.6f %.6f %.6f \n", iout, norm(sol.u[iout][:], 1), norm(sol.u[iout][:], 2), norm(sol.u[iout][:], Inf))
        end
    end #f
    println(string(" # Writing diagnostics to ASCII file:", OUTPUT_DIR, "*.dat ...  DONE ") )
end



function compute_mass!(uaux, u, mesh, metrics, ω,neqs)

    u2uaux!(uaux, u, neqs, mesh.npoin)
    mass = 0.0	
    for iel=1:mesh.nelem

        for j=1:mesh.ngl, i=1:mesh.ngl
            ip = mesh.connijk[iel,i,j]
            ωJac = ω[i]*ω[j]*metrics.Je[iel,i,j]
            ρ    = uaux[ip,1]
            mass += ρ*ωJac        
        end
     end
    return mass

end

function compute_energy!(uaux, u, mesh, metrics, ω,neqs)

    PhysConst = PhysicalConst{Float64}()	 
    u2uaux!(uaux, u, neqs, mesh.npoin)
    
    energy = 0.0	
    for iel=1:mesh.nelem
        for j=1:mesh.ngl, i=1:mesh.ngl
            ip = mesh.connijk[iel,i,j]
            ωJac = ω[i]*ω[j]*metrics.Je[iel,i,j]
            ρ    = uaux[ip,1]
            u    = uaux[ip,2]/ρ
            v    = uaux[ip,3]/ρ
            θ   = uaux[ip,4]/ρ		
            P = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
	    exner = (P/PhysConst.pref)^(1/PhysConst.cpoverR)
            T =  θ*exner  
            z = mesh.z[ip]                                        #Only valid for flow in a box
	    ke = 0.5*ρ*(u*u + v*v)
            ie = ρ*PhysConst.cv*T
            pe = ρ*PhysConst.g*z 
            te = ke + ie + pe
            energy += te* ωJac   
         end
     end
    return energy

end

function print_diagnostics(mass_ini, energy_ini, uaux, solution, mesh, metrics, ω, neqs)
    
    iout = inputs[:ndiagnostics_outputs]
    mass_final = compute_mass!(uaux, @view(solution.u[iout][:]), mesh, metrics, ω, neqs)
    energy_final = compute_energy!(uaux, @view(solution.u[iout][:]), mesh, metrics, ω, neqs)
    mass_loss = abs(mass_final-mass_ini)/mass_ini
    energy_loss = abs(energy_final-energy_ini)/energy_ini

    if mass_loss > 1.0e-11
        
        println(" # --- Mass   Loss: ", RED_FG(string(mass_loss)), RED_FG(" !!! WARNING: Large value. Maybe diffusion is on."))
        println(" # --- Energy Loss: ", energy_loss)
    else
        println(" # --- Mass   Loss: ",mass_loss)
        println(" # --- Energy Loss: ",energy_loss)
    end
    #println("linf norm of u: ", norm(u, Inf))
end
