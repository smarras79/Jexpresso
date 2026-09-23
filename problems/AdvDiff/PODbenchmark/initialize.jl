#---------------------------------------------------------------------------------
# PODbenchmark — the initial condition of the POD reference benchmark.
#
#   u(x,0) = Σ_{j=1}^3 A_j cos(2πj x/L + ϕ_j) ,   A = (1, ½, ¼) ,
#                                                 ϕ = (0, 0.7, -1.3)
#
# and the exact solution is that profile translated at speed c: u(x,t) =
# u(x - ct, 0). THE AMPLITUDES ARE THE ANSWER — the POD of this flow is known in
# closed form and is set by them alone (README, Eqs. (B2)-(B5)):
#
#   λ_{2j-1} = λ_{2j} = A_j² L/4 ,   E_j = A_j²/(2 Σ A²) ,
#
# i.e. one degenerate PAIR per wavenumber, carrying 38.10 %, 9.52 % and 2.38 %
# of the energy twice each. The phases ϕ_j do not appear: shifting a harmonic
# rotates the two modes of its pair within their own plane and changes nothing
# that is observable. They are non-zero here precisely so that a run cannot
# agree with the reference for the wrong reason — a code that quietly assumed
# the wave was a pure cosine would still get the spectrum right and the phase
# portraits wrong.
#
# Change A to change the spectrum, and the reference with it. A = (1,1,1,…)
# makes the spectrum FLAT, which is the Kolmogorov n-width statement of the
# README: POD cannot compress a travelling wave.
#---------------------------------------------------------------------------------
function initialize(SD, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" # Initialize fields for the POD benchmark (advected multi-harmonic wave) ....")

    qvars = ["u"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))

    L = Float64(inputs[:xmax] - inputs[:xmin])
    A = (1.0, 0.5, 0.25)
    ϕ = (0.0, 0.7, -1.3)

    for ip = 1:mesh.npoin
        x = mesh.x[ip] - Float64(inputs[:xmin])
        s = 0.0
        for j = 1:length(A)
            s += A[j]*cos(2π*j*x/L + ϕ[j])
        end
        q.qn[ip,1] = s
        q.qe[ip,1] = 0.0
    end

    println(" #   L = ", L, " ; A = ", A, " ; ϕ = ", ϕ)
    println(" #   exact POD: λ = ", join((string(round(a^2*L/4, digits=6), " (×2)") for a in A), ", "))
    println(" # Initialize fields for the POD benchmark .................... DONE")

    return q
end
