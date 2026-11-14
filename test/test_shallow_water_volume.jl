# :swe_conserve_volume — the volume compensation of the continuity tendency.
# Run (single rank): julia --project=. --startup-file=no test/test_shallow_water_volume.jl
using Test, MPI
MPI.Initialized() || MPI.Init()

module SWVolumeTest
    using MPI
    get_mpi_comm() = MPI.COMM_WORLD
    struct CPU end
    struct ContGal end
    include(joinpath(@__DIR__, "..", "src", "kernel", "physics", "shallow_water_volume.jl"))
end
const SWV = SWVolumeTest

volume_params(M, owner) = (M=M, mesh=(npoin=length(M), gip2owner=owner))
compensate!(RHS, p) = SWV.conserve_shallow_water_volume!(RHS, p)
alloc_of(RHS, p)    = @allocated compensate!(RHS, p)

@testset "shallow-water volume compensation" begin
    npoin = 200
    M     = 0.5 .+ rand(npoin)                 # lumped mass, all owned by rank 0
    RHS0  = randn(npoin, 3)
    RHS0[:, 1] .+= 3.0 .* M                    # a net source: d/dt ∫H dA ≠ 0
    p     = volume_params(M, zeros(Int, npoin))

    RHS = copy(RHS0)
    w = compensate!(RHS, p)
    @test w ≈ sum(RHS0[:, 1]) / sum(M)                   # the removed basin-mean rate
    @test abs(sum(RHS[:, 1])) <= 1e-12 * sum(abs, RHS0[:, 1])   # volume tendency gone
    @test RHS[:, 1] ≈ RHS0[:, 1] .- w .* M               # a uniform detrainment, nothing else
    @test RHS[:, 2:3] == RHS0[:, 2:3]                    # momenta untouched

    # nodes owned by another rank enter neither sum (they are counted there)
    owner = zeros(Int, npoin); owner[1:50] .= 1
    RHS = copy(RHS0)
    w2 = compensate!(RHS, volume_params(M, owner))
    @test w2 ≈ sum(RHS0[51:end, 1]) / sum(M[51:end])
    @test RHS[:, 1] ≈ RHS0[:, 1] .- w2 .* M

    # no allocation (after compilation)
    RHS = copy(RHS0); compensate!(RHS, p)
    @test alloc_of(RHS, p) == 0

    @test_throws ErrorException compensate!(copy(RHS0), (M=rand(4, 4), mesh=(npoin=4, gip2owner=zeros(Int, 4))))
end
