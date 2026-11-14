# Mass-conservation writer (:lmass_diagnostic): file contents and no allocation per write.
# Run (single rank): julia --project=. --startup-file=no test/test_mass_diagnostic.jl
using Test, MPI
MPI.Initialized() || MPI.Init()

module MassDiagTest
    using MPI, Printf
    get_mpi_comm() = MPI.COMM_WORLD
    struct CPU end
    struct ContGal end
    include(joinpath(@__DIR__, "..", "src", "io", "mass_diagnostic.jl"))
end
const MD = MassDiagTest

function fake_params(M, qe1)
    npoin = length(M)
    qe = zeros(npoin, 3); qe[:, 1] .= qe1
    return (M=M, mesh=(npoin=npoin, gip2owner=zeros(Int, npoin)), qp=(qe=qe,),
            backend=MD.CPU(), AD=MD.ContGal())
end
write_once!(d, integ) = MD.mass_diagnostic_affect!(d, integ)
alloc_of(d, integ)    = @allocated write_once!(d, integ)

@testset "mass diagnostic writer" begin
    dir = mktempdir()
    npoin = 300
    M  = 0.5 .+ rand(npoin)
    H0 = 150.0 .+ randn(npoin)
    u0 = vcat(H0, randn(2npoin))                   # field-major [H; Hu; Hv]
    p  = fake_params(M, 0.0)
    d  = MD.MassDiagnostic(; interval=2)
    MD.open_mass_diagnostic!(d, u0, p, Dict{Symbol,Any}(:output_dir => dir), 0.0)
    @test d.active && isfile(joinpath(dir, "mass.csv"))
    @test d.mass0 ≈ sum(M .* H0)

    # the condition fires every `interval` accepted steps
    @test [MD.mass_diagnostic_condition(d) for _ in 1:4] == [false, true, false, true]

    u1 = copy(u0); u1[1:npoin] .*= 1.001          # +0.1 % mass
    integ = (u=u1, p=p, t=7.776e6)
    write_once!(d, integ)                          # compile
    @test alloc_of(d, integ) == 0                  # a write allocates nothing
    MD.close_mass_diagnostic!(d)
    @test !d.active && d.io === nothing

    lines = readlines(joinpath(dir, "mass.csv"))
    data  = [split(l, ',') for l in lines if !startswith(l, '#') && !startswith(l, "step")]
    @test length(data) == 3                        # t0 line + two writes
    @test parse(Float64, data[1][4]) ≈ d.mass0 && parse(Float64, data[1][5]) == 0.0
    @test parse(Float64, data[2][2]) == 7.776e6
    @test parse(Float64, data[2][3]) ≈ 90.0         # t in days
    @test parse(Float64, data[2][5]) ≈ 1.0e-3 rtol=1e-6

    # perturbation form: the state holds H - He, the mass is of H
    pp = fake_params(M, 150.0)
    d2 = MD.MassDiagnostic(; lpert=true)
    MD.open_mass_diagnostic!(d2, vcat(H0 .- 150.0, zeros(2npoin)), pp, Dict{Symbol,Any}(:output_dir => mktempdir()), 0.0)
    @test d2.mass0 ≈ sum(M .* H0)
    MD.close_mass_diagnostic!(d2)
end
