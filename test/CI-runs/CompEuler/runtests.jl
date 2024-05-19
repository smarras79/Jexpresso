module CompEulerTests

using Test

@testset "theta" begin include("theta/Tests.jl") end
@testset "thetaTracers" begin include("thetaTracers/Tests.jl") end
@testset "theta_laguerre" begin include("theta_laguerre/Tests.jl") end
@testset "wave1d_lag" begin include("wave1d_lag/Tests.jl") end

end # module
