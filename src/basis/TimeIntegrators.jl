mutable struct RK3_Integrator{TFloat}
  a::Array{TFloat}
  b::Array{TFloat}
  g::Array{TFloat}
end

function buildRK3Integrator!(RK::RK3_Integrator)
  RK.a=[0.0 -5/9 -153/128]
  RK.b=[0.0 1/3 3/4]
  RK.g=[1/3 15/16 8/15]
end
