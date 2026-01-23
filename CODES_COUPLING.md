mpirun --tag-output \
    -np 2 ./alya/Alya.x : \
    -np 2 julia --project=. src/Jexpresso.jl CompEuler wave1d false \
        --gather-coupling --coupling-test-only --code-name "Jexpresso"


When you have the full Alya simulation that continues running after the gather, then you can omit --coupling-test-only and both codes will run together.

Try with --coupling-test-only now:

mpirun --tag-output \
    -np 2 ./alya/Alya.x : \
    -np 2 julia --project=. src/Jexpresso.jl CompEuler wave1d false \
        --gather-coupling --coupling-test-only --code-name "Jexpresso"