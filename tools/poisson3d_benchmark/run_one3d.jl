#=============================================================================
 run_one3d.jl — run ONE configuration of the 3D periodic Poisson deck
 (problems/Elliptic/poisson_periodic_sem_3d) through Jexpresso.run_case and
 write its solution for ParaView.

   julia --project=. -t 4 tools/poisson3d_benchmark/run_one3d.jl \
         --solver sem --ne 8 --nop 4 [--outdir vis3d] [--ordering metis|cholmod]

   --solver   sem | sem_amg | sc_direct | sc_amg   (SEM: <outdir>/…/iter_1.pvtu, field u)
              ps | fft                             (uniform (ne·N)³ grid: …/pseudospectral_laplace.vtk
                                                    or …/fft_laplace.vtk, fields u, u_exact, error)
   --ne       elements per direction (ne³ hexahedra on [0,2π]³)
   --nop      SEM order N (the Fourier solvers use the (ne·N)³ grid)
   --outdir   output root (default vis3d); the files land in
              <outdir>/Elliptic/poisson_periodic_sem_3d/output/

 The L∞ error against the exact solution is printed at the end.
=============================================================================#
using Jexpresso, Printf

function main(args)
    o = Dict{String, String}("--solver" => "sem", "--ne" => "8", "--nop" => "4",
                             "--outdir" => "vis3d", "--ordering" => "metis")
    i = 1
    while i <= length(args)
        haskey(o, args[i]) && i < length(args) || error("unknown option or missing value: $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    s = Symbol(o["--solver"]); ne = parse(Int, o["--ne"]); N = parse(Int, o["--nop"])
    s in (:sem, :sem_amg, :sc_direct, :sc_amg, :ps, :fft) || error("--solver: sem sem_amg sc_direct sc_amg ps fft")
    outdir = abspath(o["--outdir"])
    ov = Dict{Symbol, Any}(
        :nop => N, :nelx => ne, :nely => ne, :nelz => ne, :fft_N => ne * N,
        :lfft => s === :fft, :lpseudospectral => s === :ps,
        :linsolve_amg => s === :sem_amg,
        :lstatic_condensation => s in (:sc_direct, :sc_amg),
        :EL_skeleton_solver => s === :sc_amg ? "amg" : "direct",
        :sparse_ordering => o["--ordering"],
        :luse_mesh_cache => false, :lbenchmark_solve => false,
        :outformat => "vtk", :output_dir => outdir)
    Jexpresso.run_case("Elliptic", "poisson_periodic_sem_3d"; inputs = ov)
    dir = joinpath(outdir, "Elliptic", "poisson_periodic_sem_3d", "output")
    @printf("\n%s on %d³ elements, N = %d (%d unknowns): L∞ error %.3e\n", s, ne, N, (ne * N)^3,
            Jexpresso.JX_LAST_SOLVE_ERR[].linf)
    println("open in ParaView: ", s in (:ps, :fft) ?
            joinpath(dir, s === :ps ? "pseudospectral_laplace.vtk" : "fft_laplace.vtk") :
            joinpath(dir, "iter_1.pvtu"))
end

main(copy(ARGS))
