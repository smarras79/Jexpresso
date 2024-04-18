module CompEulerTests

using Test

@testset "theta" begin include("theta/thetaTests.jl") end
@testset "thetaTracers" begin include("thetaTracers/thetaTests.jl") end

end # module
