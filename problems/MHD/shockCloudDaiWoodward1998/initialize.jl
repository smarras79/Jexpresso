#---------------------------------------------------------------------------------
# Interaction of a Mach 12 MHD fast shock with dense clouds shaped as the letters "MHD":
# W. Dai & P. R. Woodward, JCP 142:331-369 (1998), p. 360 and Fig. 20. Domain [0, 2] × [0, 1.5], γ = 5/3.
#
#   pre-shock  (ρ, p, u, v, w, Bx, By, Bz) = (1, 1, 0, 0, 0, 2, 2, 2)                          x > x_s
#   post-shock (ρ, p, u, v, w, Bx, By, Bz) = (3.91109, 253.724, 13.8388, −0.0500186, −0.0500186,
#                                             2, 7.84323, 7.84323)                              x < x_s
#   letters    ρ = 50, the pre-shock p, v and B (pressure equilibrium, at rest)
#
# B as printed by DW is Gaussian (B²/8π); it enters here divided by √(4π) (Heaviside-Lorentz, ½|B|²),
# the conversion Tóth (2000, §6.5) applies to DW's single-cloud test. With it the two states satisfy the
# Rankine-Hugoniot conditions to 6 digits with shock speed s = ρ₂u₂/(ρ₂ − ρ₁) = 18.59 (README).
#
# Neither the shock position x_s nor the letters' geometry is printed. Both are read off Fig. 20:
#   x_s = 0.2      the incident shock is at x = 0.664 at t = 0.02507 and 18.59 × 0.02507 = 0.466
#   letters        Times Bold, cap height 0.395, from x = 0.243 to 1.680 (mhd_letters.jl)
#
# :sc_smooth (default 0, DW's sharp data) spreads the letters' edge and the shock over that width with a
# C² ramp of the signed distance (letters) or of x − x_s (shock), blending the primitive variables.
#---------------------------------------------------------------------------------

# Even-odd point-in-polygon test over all the letters' contours, and the distance to the nearest edge.
function sc_letters_signed_distance(x, y)
    inside = false
    d2 = Inf
    for letter in SC_LETTERS, c in letter
        n = size(c, 2)
        j = n
        for i = 1:n
            xi, yi, xj, yj = c[1,i], c[2,i], c[1,j], c[2,j]
            if (yi > y) != (yj > y) && x < (xj - xi)*(y - yi)/(yj - yi) + xi
                inside = !inside
            end
            ex, ey = xj - xi, yj - yi
            τ  = clamp(((x - xi)*ex + (y - yi)*ey)/(ex*ex + ey*ey), 0.0, 1.0)
            dx, dy = x - xi - τ*ex, y - yi - τ*ey
            d2 = min(d2, dx*dx + dy*dy)
            j = i
        end
    end
    return inside ? sqrt(d2) : -sqrt(d2)       # > 0 inside a letter
end

# Weight 0 → 1 of the second state across a jump at signed distance d: a step for w ≤ 0, else a C² ramp over |d| < w/2.
@inline function sc_ramp(d, w)
    w <= 0.0 && return d > 0.0 ? 1.0 : 0.0
    ξ = clamp(d/w + 0.5, 0.0, 1.0)
    return ξ*ξ*ξ*(10.0 + ξ*(-15.0 + 6.0*ξ))
end

function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank == 0 && @info " Initialize fields for 2D ideal GLM-MHD (Dai-Woodward shock / MHD clouds) ........... "

    qvars    = ["ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ"]
    qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T", "pmag", "Mach", "β", "Mach_f", "log10ρ"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    inputs[:backend] != CPU()        && error(" problems/MHD/shockCloudDaiWoodward1998: only the CPU backend is supported.")
    inputs[:SOL_VARS_TYPE] != TOTAL() && error(" problems/MHD/shockCloudDaiWoodward1998: only SOL_VARS_TYPE = TOTAL() is supported.")

    xs = inputs[:sc_shock_x]::Float64
    w  = inputs[:sc_smooth]::Float64
    pre, post = SC_PRE, SC_POST
    ncloud = 0
    for ip = 1:mesh.npoin
        x, y = mesh.coords[1,ip], mesh.coords[2,ip]
        a = sc_ramp(sc_letters_signed_distance(x, y), w)    # 1 inside the letters
        b = sc_ramp(xs - x, w)                              # 1 behind the shock
        # pre-shock gas with the letters' density, then the post-shock state blended in (primitives)
        prim = ntuple(k -> k == 1 ? (1.0 - a)*pre[1] + a*SC_RHO_CLOUD : pre[k], 8)
        prim = ntuple(k -> (1.0 - b)*prim[k] + b*post[k], 8)
        s = sc_conserved(prim)
        for ieq = 1:9
            q.qn[ip,ieq] = s[ieq]
            q.qe[ip,ieq] = s[ieq]
        end
        q.qn[ip,end] = prim[2]
        q.qe[ip,end] = prim[2]
        ncloud += a > 0.5
    end
    # GLM cleaning speed: the post-shock gas is the fastest, |v| + c_f = 24.4
    c_h_mhd[] = max(sc_wave_speed(pre...), sc_wave_speed(post...))

    if rank == 0
        (abs(mesh.xmin) > 1e-8 || abs(mesh.xmax - 2.0) > 1e-8 || abs(mesh.ymin) > 1e-8 || abs(mesh.ymax - 1.5) > 1e-8) &&
            @warn " problems/MHD/shockCloudDaiWoodward1998: the letters are placed in [0, 2] × [0, 1.5], but this mesh spans [$(mesh.xmin), $(mesh.xmax)] × [$(mesh.ymin), $(mesh.ymax)]."
        @info " shock at x = $(xs), edge smoothing width $(w), $(ncloud) nodes inside the letters (rank 0)"
        @info " GLM divergence-cleaning speed c_h = $(c_h_mhd[])"
        @info " Initialize fields for 2D ideal GLM-MHD (Dai-Woodward shock / MHD clouds) ........... DONE"
    end
    return q
end
