module CompEulerTests

using Test

@testset "theta" begin include("theta/thetaTests.jl") end
@testset "thetaTracers" begin include("thetaTracers/thetaTests.jl") end
@testset "wave1d_lag" begin include("wave1d_lag/Tests.jl") end

end # module
