# ASCII
function solution_norms(sol::ODESolution, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII; )
    
    println(string(" # Writing diagnostics to ASCII file:", OUTPUT_DIR, "*.dat ...  ") )

    fname = @sprintf "diagnostics.dat"
    open(string(OUTPUT_DIR, "/", fname), "w") do f
        @printf(f, "time |  norm₁ | norm₂ | norm∞ \n")
        for iout = 1: inputs[:ndiagnostics_outputs]
            @printf(f, " %.4f %.6f %.6f %.6f \n", sol.t, norm(sol.u[iout][:], 1), norm(sol.u[iout][:], 1), norm(sol.u[iout][:], Inf))
        end
    end #f     
    println(string(" # Writing diagnostics to ASCII file:", OUTPUT_DIR, "*.dat ...  DONE ") )
end

