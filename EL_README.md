# Steps to run the element learning case

julia> push!(empty!(ARGS), "Elliptic", "ElementLearning")

1)
1st run:
lEL_Train      => true,
:gmsh_filename => "./meshes/gmsh_grids/square_dirichletT_1x1.msh",

julia> include("./src/Jexpresso.jl")

2)
cd EL_Jexpresso
Set nepochs in train.py
./training_script.sh

3)
2nd run
lEL_Train      => false,
:gmsh_filename => "./meshes/gmsh_grids/square_dirichletT_15x15.msh",
