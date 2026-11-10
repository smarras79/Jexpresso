#---------------------------------------------------------------------------------
# Shock / dense cloud interactions of W. Dai & P. R. Woodward, JCP 142:331-369 (1998), γ = 5/3.
# The states are in user_flux.jl (SC_STATES); B as printed by DW is Gaussian (B²/8π) and enters divided
# by √(4π) (Heaviside-Lorentz, ½|B|²), the conversion Tóth (2000, §6.5) applies to DW's single-cloud test.
# With it the states of both tests satisfy the Rankine-Hugoniot conditions to 5-6 digits (README).
#
# :sc_problem => "letters"  (p. 360, Fig. 20) Mach 12 shock moving right at s = 18.59 into gas at rest;
#   [0, 2] × [0, 1.5]; post-shock state for x < x_s; ahead of it the letters "MHD" with ρ = 50 and the
#   pre-shock p, v, B. Not printed by DW, read off Fig. 20:
#     x_s = 0.2     the incident shock is at x = 0.664 at t = 0.02507, and 18.59 × 0.02507 = 0.466
#     letters       Times Bold, cap height 0.395, from x = 0.243 to 1.680 (SC_LETTERS below)
#
# :sc_problem => "cloud"    (p. 365, Figs. 23-25) Mach 15 shock; the shocked gas is at rest and the gas
#   ahead, cloud included, flows in from the right at u = −15.27 (shock speed +5.18 in this frame);
#   [0, 1.5] × [0, 1]; a round cloud with ρ = 10 and the pre-shock p, v, B. Its centre and radius are not
#   printed for this test; they are those DW print for their Mach 10 cloud (Fig. 18), with the shock where
#   Tóth (2000, §6.5) puts it for that test:
#     x_s = 0.6,  cloud centre (0.8, 0.5),  radius 0.15
#
# Both clouds are defined analytically: a node's density is set from its coordinates (point-in-polygon
# for the letters, distance to the centre for the round cloud), on any mesh.
#
# :sc_smooth spreads the shock and the cloud edges over that width with a C² ramp of x − x_s (shock) or of
# the signed distance to the edge (clouds), blending the primitive variables; 0 gives DW's sharp data.
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# The letters "MHD" of Fig. 20 as polygons in the domain's coordinates (x in row 1, y in row 2).
# The letters of Fig. 20 are Times Bold. These vertices are the outlines of M, H, D in Liberation Serif
# Bold (Times-metric), quadratic arcs flattened into 6 segments each, scaled to the cap height measured on
# the figure's first panel (t = 0.02507; base y = 0.5746, cap line y = 0.9697) and placed where that panel
# has them: H from x = 0.8142, D from x = 1.2839, M ending at x = 0.7850 (its left side is already swept
# by the shock in that panel). Overlaid on the scan they match DW's H and D to a pixel.
# One tuple of contours per letter; D's second contour is its counter. Inside = odd number of crossings.
#---------------------------------------------------------------------------------
if !@isdefined(SC_LETTERS)
    const SC_LETTERS = (
        # M
        (
            [0.49273 0.47652 0.33068 0.33068 0.38371 0.38371 0.24317 0.24317 0.29385 0.29385 0.24317 0.24317 0.39844 0.51129 0.62649 0.78500 0.78500 0.73432 0.73432 0.78500 0.78500 0.58848 0.58848 0.64151 0.64151;
             0.57460 0.57460 0.90842 0.60406 0.59611 0.57460 0.57460 0.59611 0.60406 0.94053 0.94819 0.96970 0.96970 0.70984 0.96970 0.96970 0.94819 0.94053 0.60406 0.59611 0.57460 0.57460 0.59611 0.60406 0.90842],
        ),
        # H
        (
            [0.81420 0.81420 0.86488 0.86488 0.81420 0.81420 1.00836 1.00836 0.95769 0.95769 1.11914 1.11914 1.06847 1.06847 1.26322 1.26322 1.21225 1.21225 1.26322 1.26322 1.06847 1.06847 1.11914 1.11914 0.95769 0.95769 1.00836 1.00836;
             0.57460 0.59640 0.60406 0.94024 0.94819 0.96970 0.96970 0.94819 0.94024 0.79410 0.79410 0.94024 0.94819 0.96970 0.96970 0.94819 0.94024 0.60406 0.59640 0.57460 0.57460 0.59640 0.60406 0.76169 0.76169 0.60406 0.59640 0.57460],
        ),
        # D
        (
            [1.58177 1.58090 1.57830 1.57396 1.56789 1.56008 1.55054 1.53904 1.52537 1.50951 1.49148 1.47128 1.44889 1.42709 1.42709 1.43637 1.44536 1.45405 1.46245 1.47055 1.47836 1.49068 1.50211 1.51264 1.52229 1.53104 1.53890 1.54599 1.55241 1.55816 1.56326 1.56769 1.57146 1.57461 1.57719 1.57919 1.58063 1.58149;
             0.77200 0.79922 0.82391 0.84607 0.86570 0.88280 0.89737 0.90957 0.91955 0.92731 0.93285 0.93618 0.93729 0.93729 0.60819 0.60747 0.60688 0.60642 0.60609 0.60590 0.60583 0.60629 0.60766 0.60996 0.61316 0.61729 0.62233 0.62835 0.63541 0.64351 0.65264 0.66282 0.67404 0.68654 0.70057 0.71613 0.73323 0.75185],
            [1.47158 1.50549 1.53632 1.56406 1.58871 1.61028 1.62877 1.64429 1.65700 1.66688 1.67394 1.67818 1.67959 1.67824 1.67419 1.66744 1.65798 1.64583 1.63098 1.61339 1.59303 1.56991 1.54403 1.51537 1.48395 1.35167 1.28390 1.28390 1.33458 1.33458 1.28390 1.28390;
             0.96970 0.96836 0.96435 0.95766 0.94829 0.93625 0.92153 0.90406 0.88377 0.86065 0.83471 0.80595 0.77436 0.74192 0.71237 0.68571 0.66194 0.64106 0.62307 0.60790 0.59549 0.58583 0.57894 0.57480 0.57342 0.57460 0.57460 0.59611 0.60406 0.94053 0.94819 0.96970],
        ),
    )
end

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

# Signed distance to the round cloud of the "cloud" test (> 0 inside).
@inline sc_disc_signed_distance(x, y) = 0.15 - sqrt((x - 0.8)^2 + (y - 0.5)^2)

function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    pb   = inputs[:sc_problem]::String
    haskey(SC_STATES, pb) || error(" problems/MHD/shockCloudDaiWoodward1998: :sc_problem = \"$(pb)\", expected \"letters\" or \"cloud\".")
    rank == 0 && @info " Initialize fields for 2D ideal GLM-MHD (Dai-Woodward shock / $(pb == "letters" ? "MHD letters" : "round cloud")) ........... "

    qvars    = ["ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ"]
    qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T", "pmag", "Mach", "β", "Mach_f", "log10ρ"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    inputs[:backend] != CPU()        && error(" problems/MHD/shockCloudDaiWoodward1998: only the CPU backend is supported.")
    inputs[:SOL_VARS_TYPE] != TOTAL() && error(" problems/MHD/shockCloudDaiWoodward1998: only SOL_VARS_TYPE = TOTAL() is supported.")

    xs = inputs[:sc_shock_x]::Float64
    w  = inputs[:sc_smooth]::Float64
    st = SC_STATES[pb]
    pre, post = st.pre, st.post
    dist = pb == "letters" ? sc_letters_signed_distance : sc_disc_signed_distance
    ncloud = 0
    for ip = 1:mesh.npoin
        x, y = mesh.coords[1,ip], mesh.coords[2,ip]
        a = sc_ramp(dist(x, y), w)       # 1 inside the cloud
        b = sc_ramp(xs - x, w)           # 1 behind the shock
        # pre-shock gas with the cloud's density, then the post-shock state blended in (primitives)
        prim = ntuple(k -> k == 1 ? (1.0 - a)*pre[1] + a*st.ρcloud : pre[k], 8)
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
    # GLM cleaning speed: the larger |v| + c_f of the two states (24.4 letters, 16.7 cloud); the cloud,
    # denser at the same p and B, is slower than the gas around it
    c_h_mhd[] = max(sc_wave_speed(pre...), sc_wave_speed(post...))
    # characteristic far-field boundaries (user_bc.jl): needs c_h
    h = sqrt((mesh.xmax - mesh.xmin)*(mesh.ymax - mesh.ymin)/mesh.gnelem)   # element size (square elements, global count)
    s = sc_build_far_field!(st, xs, 2.0*h + w)

    if rank == 0
        Lx, Ly = pb == "letters" ? (2.0, 1.5) : (1.5, 1.0)
        (abs(mesh.xmin) > 1e-8 || abs(mesh.xmax - Lx) > 1e-8 || abs(mesh.ymin) > 1e-8 || abs(mesh.ymax - Ly) > 1e-8) &&
            @warn " problems/MHD/shockCloudDaiWoodward1998: the \"$(pb)\" test is set in [0, $(Lx)] × [0, $(Ly)], but this mesh spans [$(mesh.xmin), $(mesh.xmax)] × [$(mesh.ymin), $(mesh.ymax)]."
        @info " shock at x = $(xs) moving at s = $(s), edge smoothing width $(w), $(ncloud) nodes inside the cloud (rank 0)"
        @info " GLM divergence-cleaning speed c_h = $(c_h_mhd[])"
        @info " Initialize fields for 2D ideal GLM-MHD (Dai-Woodward shock) ........... DONE"
    end
    return q
end
