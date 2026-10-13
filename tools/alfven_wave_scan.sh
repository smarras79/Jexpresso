#!/usr/bin/env bash
#
# Convergence history of the smooth Alfvén wave (Babbar & Ranocha 2026, §5.5.1
# and Fig. 13): every polynomial order on every mesh, in ONE Julia session (the
# case is compiled once), which problems/MHD/alfvenWaveBabbarRanocha2026/
# user_plot.jl turns into output/convergence-it<n>.png and a table on stdout.
#
#   tools/alfven_wave_scan.sh                       # N = 2 3 on 16² ... 128², DSGS and Galerkin
#   AW_NOPS="3 4" AW_NELX="8 16 32" tools/alfven_wave_scan.sh
#   AW_VISC=dsgs tools/alfven_wave_scan.sh          # only the DynSGS runs
#   JULIA=/path/to/julia tools/alfven_wave_scan.sh
#
# The defaults are the paper's meshes, 16², 32², 64², 128² elements, and its
# orders except N = 1: Jexpresso cannot run N = 1 on a periodic 2D mesh yet
# (add_high_order_nodes_edges! returns before filling poin_in_edge when
# nop < 2, and the periodic setup of mod_mesh_read_gmsh! then reads node 0;
# BoundsError at src/kernel/mesh/mesh.jl:3456). A run that fails is reported
# and the sweep goes on with the next one.
#
# Both the DynSGS run and the plain Galerkin one are made by default
# (AW_VISC="dsgs none"): on a smooth solution the residual viscosity must not
# cost the scheme its order.
#
# THE ERROR STORE IS CLEARED FIRST, since every run redraws the figure from all
# the errors stored for the same final time; AW_KEEP=1 accumulates instead, to
# finish a partial sweep.
set -u
cd "$(dirname "$0")/.."

NOPS=${AW_NOPS:-"2 3"}
NELX=${AW_NELX:-"16 32 64 128"}
VISCS=${AW_VISC:-"dsgs none"}
JULIA=${JULIA:-julia}

[ "${AW_KEEP:-0}" = "1" ] || rm -rf problems/MHD/alfvenWaveBabbarRanocha2026/errors
AW_NELX="$NELX" tools/alfven_wave_mesh.sh

GKSwstype=${GKSwstype:-100} "$JULIA" --project=. -e '
using Jexpresso
nops, nelxs, viscs = (split(a) for a in ARGS)
failed = String[]
for v in viscs, n in nops, m in nelxs
    ENV["JEXPRESSO_AW_VISC"] = v; ENV["JEXPRESSO_AW_NOP"] = n; ENV["JEXPRESSO_AW_NELX"] = m
    println("=== smooth Alfvén wave: N = ", n, ", ", m, "x", m, " elements, visc ", v)
    t0 = time()
    try
        Jexpresso.run_case("MHD", "alfvenWaveBabbarRanocha2026")
        println("=== N = ", n, ", ", m, "x", m, ", ", v, ": ", round(time() - t0, digits = 1), " s")
    catch err
        push!(failed, string("N = ", n, ", ", m, "x", m, ", ", v))
        println("=== N = ", n, ", ", m, "x", m, ", ", v, ": FAILED, ", first(sprint(showerror, err), 300))
    end
end
println("=== done; the figure and the table of the last run hold the whole sweep",
        isempty(failed) ? "" : string(". FAILED: ", join(failed, "; ")))
' "$NOPS" "$NELX" "$VISCS"
