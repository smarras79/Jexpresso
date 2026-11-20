#-------------------------------------------------------------------------------------------------
#   ∂u/∂t = ν ∂²u/∂x²  on [0, L],   ∂u/∂n = 0 at both ends (insulated)
#
# In the weak form the boundary term is ∮ ψ ν ∂u/∂n dS. Leaving it out, which is what the kernel
# does when no flux is supplied, IS the zero-flux condition. So an insulated wall needs no
# user_bc_neumann!, no :bdy_fluxes and nothing in user_bc_dirichlet!. Every other condition
# starts from here.
#
# Exact solution (two cosine modes, each with zero slope at both ends):
#   u(x, t) = 1 + A1 exp(-ν (π/L)² t) cos(π x/L) + A2 exp(-ν (2π/L)² t) cos(2π x/L)
# The mean of u is conserved: nothing enters or leaves.
#-------------------------------------------------------------------------------------------------
const I1_NU = Ref(0.1)
const I1_L  = Ref(1.0)
const I1_A1 = 0.5
const I1_A2 = 0.25

function i1_exact(x, t)
    k1 = π / I1_L[]; k2 = 2π / I1_L[]
    return 1 + I1_A1 * exp(-I1_NU[] * k1^2 * t) * cos(k1 * x) + I1_A2 * exp(-I1_NU[] * k2^2 * t) * cos(k2 * x)
end

function initialize(SD::NSD_1D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)
    I1_NU[] = Float64(inputs[:μ][1])
    I1_L[]  = Float64(inputs[:xmax] - inputs[:xmin])
    qvars = ["q"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs = length(qvars))
    for ip = 1:mesh.npoin
        q.qn[ip, 1] = i1_exact(mesh.coords[1, ip] - inputs[:xmin], 0.0)
        q.qe[ip, 1] = 0.0
    end
    return q
end
