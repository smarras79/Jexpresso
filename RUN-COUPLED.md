
#
# NEW WITH POSSIBLY FIXED CHART COMMUNICATION:
#
mpirun --tag-output -np 2 ./alya/Alya_enhanced.x : \
       -np 2 julia --project=. Jexpresso-mini-coupled.jl

#
# OLD POSSIBLY INCORRECT
#
#mpirun --tag-output -np 2 ./alya/Alya.x : -np 2 julia --project=. Jexpresso-mini-coupled.jl false --gather-coupling --coupling-test-only --code-name "Jexpresso"