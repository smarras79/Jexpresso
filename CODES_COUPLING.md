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

Scenario	Command	Use Case
Quick coupling test	coupling_test_minimal.jl	Testing with Alya unit tests (recommended)
Test after module loading	--coupling-test-only	Verify Jexpresso modules load in MPI environment
Full simulation	--gather-coupling (no test flag)	Production coupling with full Alya simulation