# Advection (speed 1) + diffusion (nu) of a Gaussian on the infinite line:
# u(x,t) = s/sqrt(s^2+2 nu t) exp(-(x-x0-t)^2/(2(s^2+2 nu t))).
# q = u(x,0); qe = u(x,tend), kept only for comparison (PERT, unused by the RHS).
function initialize(SD, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)
    qvars = ("q")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))
    # x0 such that q(x = 0) = exp(-x0^2/18) ~ 3e-18: the inflow value q = 0 then
    # matches the exact solution to machine precision (at x0 = 20 the mismatch,
    # 2e-10, set the error floor)
    x0 = 27.0; s = 3.0
    nu = inputs[:lvisc] ? inputs[:μ][1] : 0.0
    tf = inputs[:tend]
    for ip=1:mesh.npoin
        x = mesh.coords[1,ip]
        q.qn[ip,1] = exp(-(x-x0)^2/(2*s^2))
        s2 = s^2 + 2*nu*tf
        q.qe[ip,1] = s/sqrt(s2)*exp(-(x-x0-tf)^2/(2*s2))
    end
    return q
end
