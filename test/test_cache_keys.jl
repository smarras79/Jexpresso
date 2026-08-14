#---------------------------------------------------------------------------------
# test/test_cache_keys.jl — the mesh / SEM preprocess cache is keyed per CASE.
#
#   julia --project=. test/test_cache_keys.jl
#
# The bug this pins down: running `ShallowWater/SoliWaveIsland_amr` right after
# `ShallowWater/SoliWaveIsland` reused the non-AMR run's cached mesh and metric
# terms — same .msh, same :nop — and the mismatch surfaced much later as a
# `DimensionMismatch: array could not be broadcast to match destination` from
# params_setup, with nothing in the message pointing at the cache.
#
# Nothing here touches the filesystem or MPI collectives: the helpers under
# test are pure functions of the inputs dict (plus the local rank, which is 0
# in this serial test).
#---------------------------------------------------------------------------------
using Test
using Jexpresso
using MPI

using Jexpresso: _case_tag, _cache_dir, _mesh_cache_path, _preprocess_cache_path,
                 _adaptive_mesh_run, _cache_fingerprint, _cache_fingerprint_matches,
                 _struct_schema_signature

MPI.Initialized() || MPI.Init()

# Two decks that are as similar as two decks get: same equations, same grid
# file, same polynomial order. Only the case name and the AMR switches differ.
const NONAMR = Dict{Symbol,Any}(
    :nop                => 4,
    :nsd                => 2,
    :backend            => Jexpresso.CPU(),
    :lread_gmsh         => true,
    :gmsh_filename      => "./meshes/gmsh_grids/SoliWaveIsland.msh",
    :lamr               => false,
    :ladapt             => false,
    :lpreadapt          => false,
    :linitial_refine    => false,
    :lrestart_amr       => false,
    :_parsed_equations  => "ShallowWater",
    :_parsed_case_name  => "SoliWaveIsland",
    :_case_dir          => "/repo/problems/ShallowWater/SoliWaveIsland",
)

const AMR = merge(NONAMR, Dict{Symbol,Any}(
    :_parsed_case_name  => "SoliWaveIsland_amr",
    :_case_dir          => "/repo/problems/ShallowWater/SoliWaveIsland_amr",
    :lamr               => true,
    :ladapt             => true,
    :lpreadapt          => true,
    :preadapt_max_level => 1,
    :amr_max_level      => 2,
))

@testset "cache identity" begin

    @testset "case tag" begin
        @test _case_tag(NONAMR) == "ShallowWater_SoliWaveIsland"
        @test _case_tag(AMR)    == "ShallowWater_SoliWaveIsland_amr"
        # Driven outside run.jl: fall back to the case directory.
        @test _case_tag(Dict{Symbol,Any}(:_case_dir => "/repo/problems/ShallowWater/TC2")) ==
              "ShallowWater_TC2"
        # No identity at all — callers must treat this as "do not cache".
        @test _case_tag(Dict{Symbol,Any}()) == ""
    end

    @testset "paths are per case" begin
        mesh_path = _mesh_cache_path(NONAMR, 1)
        sem_path  = _preprocess_cache_path(NONAMR, 4, 4, 1)
        @test occursin(joinpath("SoliWaveIsland", ".jexpresso_cache"), mesh_path)
        @test occursin("ShallowWater_SoliWaveIsland", basename(mesh_path))
        @test occursin("ShallowWater_SoliWaveIsland", basename(sem_path))

        # The same deck without :_case_dir must still not land in the shared
        # gmsh directory, where every case reading that .msh would collide.
        loose_a = delete!(copy(NONAMR), :_case_dir)
        loose_b = delete!(copy(merge(AMR, Dict{Symbol,Any}(:lamr => false, :ladapt => false,
                                                           :lpreadapt => false))), :_case_dir)
        @test !endswith(dirname(_mesh_cache_path(loose_a, 1)), "gmsh_grids")
        @test _mesh_cache_path(loose_a, 1) != _mesh_cache_path(loose_b, 1)

        # Unidentifiable case → no caching at all.
        @test _mesh_cache_path(Dict{Symbol,Any}(:nop => 4,
                                                :gmsh_filename => "./m.msh"), 1) == ""
        # Explicit opt-out still honoured.
        @test _mesh_cache_path(merge(NONAMR, Dict{Symbol,Any}(:luse_mesh_cache => false)), 1) == ""
    end

    @testset "adaptive runs never use a cache" begin
        @test !_adaptive_mesh_run(NONAMR)
        @test _adaptive_mesh_run(AMR)
        for k in (:lamr, :ladapt, :lpreadapt, :linitial_refine, :lrestart_amr)
            @test _adaptive_mesh_run(merge(NONAMR, Dict{Symbol,Any}(k => true)))
        end
        @test _mesh_cache_path(AMR, 1) == ""
        @test _preprocess_cache_path(AMR, 4, 4, 1) == ""
    end

    @testset "fingerprint" begin
        fp = _cache_fingerprint(NONAMR, 1)
        @test _cache_fingerprint_matches(fp, NONAMR, 1)
        # The non-AMR cache must not validate for the AMR case…
        @test !_cache_fingerprint_matches(fp, AMR, 1)
        # …nor for a different rank count.
        @test !_cache_fingerprint_matches(fp, NONAMR, 2)
        @test fp["__case_tag__"] == "ShallowWater_SoliWaveIsland"

        # Every input that changes the mesh or the metric terms invalidates.
        for k in (:nop, :lpreadapt, :preadapt_max_level, :amr_max_level, :lrestart_amr,
                  :lread_gmsh, :nelx, :nely, :nelz,
                  :xmin, :xmax, :ymin, :ymax, :zmin, :zmax,
                  :lspherical_shell, :sphere_metrics, :nop_laguerre)
            changed = _cache_fingerprint(merge(NONAMR, Dict{Symbol,Any}(k => 7)), 1)
            @test !_cache_fingerprint_matches(changed, NONAMR, 1)
        end

        # Adding a field to St_mesh / St_metrics must invalidate every cache on
        # disk without anyone bumping a version by hand: the signature of the
        # cached structs' field sets is part of the fingerprint.
        @test haskey(fp, "__struct_schema__")
        @test fp["__struct_schema__"] == _struct_schema_signature()
        stale = copy(fp)
        stale["__struct_schema__"] = "not-the-current-field-set"
        @test !_cache_fingerprint_matches(stale, NONAMR, 1)
    end
end
