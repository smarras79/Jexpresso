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
