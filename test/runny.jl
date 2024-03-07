using Test
using Jexpresso

function run_example(problem_name::String, case_name::String)
    ENV["JEXPRESSO_HOME"] = joinpath(@__DIR__, "..") 
    example_dir = joinpath(ENV["JEXPRESSO_HOME"], "test","reference", "problems", "equations", problem_name, case_name)
    @testset "$problem_name - $case_name" begin
        cd(example_dir)
        empty!(ARGS) # Clear ARGS to ensure clean state
        push!(ARGS, problem_name, case_name)
        include(joinpath(ENV["JEXPRESSO_HOME"], "src", "Jexpresso.jl"))
    end
end

# Run test sets for each example
@testset "JEXPRESSO Examples" begin
    # List of (problem_name, case_name) tuples
    examples = [
        #("CompEuler", "2d"),  ##errore in q_define sembrebbe initialize ln11
        #("CompEuler", "dc"),   # Initial Mass  :   1.458934183926933e8 and then stops
        #("CompEuler", "dc-mount"), #LoadError: MethodError: no method matching filter!
        #("CompEuler", "HSmount"), #LoadError: MethodError: no method matching DSS_rhs!(::SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}}, true}, ::SubArray{Float64, 4, Array{Float64, 4}, NTuple{4, Base.Slice{Base.OneTo{Int64}}}, true}, ::Main.Jexpresso.St_mesh{Int64, Float64}, ::Int64, ::Int64, ::Int64, ::Main.Jexpresso.NSD_2D)
        #("CompEuler", "HSmount_Lag_working"), # LoadError: BoundsError: attempt to access 47519-element Vector{Float64} at index [47520]
        #("CompEuler", "HSmount_standard"), # LoadError: MethodError: no method matching DSS_rhs!
        #("CompEuler", "mountain"), #LoadError: MethodError: no method matching define_q(::Main.Jexpresso.NSD_2D, ::Int64, ::Int64, ::Int64, ::Type{Float64}; neqs::Int64)
        #("CompEuler", "NHSmount_Lag_working"), # LoadError: BoundsError: attempt to access 79933-element Vector{Float64} at index [79934]
        #("CompEuler", "NHSmount_standard"), #LoadError: MethodError: no method matching DSS_rhs!(::SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}}, true}, ::SubArray{Float64, 4, Array{Float64, 4}, NTuple{4, Base.Slice{Base.OneTo{Int64}}}, true}, ::Main.Jexpresso.St_mesh{Int64, Float64}, ::Int64, ::Int64, ::Int64, ::Main.Jexpresso.NSD_2D)
        #("CompEuler", "NLHSmount"), # LoadError: BoundsError: attempt to access 87143-element Vector{Float64} at index [87144]
        #("CompEuler", "nozzleanderson"),
        #("CompEuler", "ScharMount"), # LoadError: MethodError: no method matching DSS_rhs!(::SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}}, true}, ::SubArray{Float64, 4, Array{Float64, 4}, NTuple{4, Base.Slice{Base.OneTo{Int64}}}, true}, ::Main.Jexpresso.St_mesh{Int64, Float64}, ::Int64, ::Int64, ::Int64, ::Main.Jexpresso.NSD_2D)
        #("CompEuler", "ScharMount_Lag"), # LoadError: BoundsError: attempt to access 16613-element Vector{Float64} at index [16614]
        #("CompEuler", "theta_pert"), LoadError: MethodError: no method matching user_flux!(::SubArray{Float64, 1, Array{Float64, 3}, Tuple{Int64, Int64, Base.Slice{Base.OneTo{Int64}}}, true}, ::SubArray{Float64, 1, Array{Float64, 3}, Tuple{Int64, Int64, Base.Slice{Base.OneTo{Int64}}}, true}, ::Main.Jexpresso.NSD_2D, ::SubArray{Float64, 1, Matrix{Float64}, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}}, true}, ::SubArray{Float64, 1, Matrix{Float64}, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}}, true}, ::Main.Jexpresso.St_mesh{Int64, Float64}, ::Main.Jexpresso.CL, ::Main.Jexpresso.PERT; neqs::Int64, ip::Int64)
        #("CompEuler", "thetaNC"),LoadError: MethodError: no method matching define_q(::Main.Jexpresso.NSD_2D, ::Int64, ::Int64, ::Int64, ::Type{Float64}; neqs::Int64)
        #("CompEuler", "thetaTracers"),
        #("CompEuler", "wave1d"),
        #("CompEuler", "wave1d_lag"),
        #("AdvDiff", "2d"),  LoadError: BoundsError: attempt to access 5×5×2 Array{Float64, 3} at index [1, 1, 3]
        #("AdvDiff", "2d_Laguerre"),
        #("AdvDiff", "2D_Wave_Train"),
        #("AdvDiff", "case1"),
        #("AdvDiff", "circle"),LoadError: periodicity requested but boundaries cannot match any vectors
        #("AdvDiff", "fd1d"),
        #("AdvDiff", "Simple_Wave"),
        #("AdvDiff", "Wave_Train_Overlapping_Plot"),
        #("AdvDiff", "Wave_Train"),
        #("Helmholtz", "case1"),
        #("Elliptic", "case1"),
    ]

    for (problem_name, case_name) in examples
        run_example(problem_name, case_name)
    end
end