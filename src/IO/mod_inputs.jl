using Crayons.Box
using Revise

export mod_inputs_user_inputs
export mod_inputs_print_welcome

include("./user_inputs.jl")

function mod_inputs_user_inputs()

    error_flag::Int8 = 0
    
    inputs = user_inputs() # user_inputs is a Dict

    print(GREEN_FG(" # User inputs: ........................ \n"))

    #
    # Check that necessary inputs exist in the Dict inside .../IO/user_inputs.jl
    #    
    if (!haskey(inputs, :equation_set))
        @error " :equation_set is missing in .../IO/user_inputs.jl"
        error_flag = 1
    else
        println( " # :equation_set =>     ", inputs[:equation_set])
    end
    
    if (!haskey(inputs, :problem))
        @error " :problem is missing in .../IO/user_inputs.jl"
        error_flag = 1
    else
        println( " # :problem      =>     ", inputs[:problem])        
    end
    
    if (haskey(inputs, :lread_gmsh) && inputs[:lread_gmsh] == true) #Non mandatory input
        println( " # lread_gmsh    =>     ", inputs[:lread_gmsh])
        if (!haskey(inputs, :gmsh_filename))
            @error " :gmsh_filename is missing in .../IO/user_inputs.jl
                     SOLUTION: add    :gmsh_filename => \"FILENAME.gmsh\" to ...IO/user_inputs.jl"
            error_flag = 1
        end
    end
    
    if (!haskey(inputs, :nop))
        @warn " :nop is missing in .../IO/user_inputs.jl
        The default value nop=4 will be used."
        inputs[:nop] = Int8(4)
        println( " # :nop          =>     ", inputs[:nop])
    else
        println( " # :nop          =>     ", inputs[:nop])
    end
    
    
    #Grid entries:
    if(!haskey(inputs, :lread_gmsh) || inputs[:lread_gmsh] == false)

        if (!haskey(inputs, :nsd))
            @error " :nsd is missing in .../IO/user_inputs.jl"
            error_flag = 1
        else
            println( " # :nsd          =>     ", inputs[:nsd])
        end
        
        if (!haskey(inputs, :npx))       
            @warn " :npx is missing in .../IO/user_inputs.jl
            The default value npx=2 will be used."
            inputs[:npx] = Int8(2)
            println( " # :npx          =>     ", inputs[:npx])
        else
            println( " # :npx          =>     ", inputs[:npx])
        end
        
        if(!haskey(inputs, :xmin))
            @warn " :xmin is missing in .../IO/user_inputs.jl
            The default value :xmin=-1.0 will be used."
            inputs[:xmin] = Float64(-1.0)
        end
        println( " # :xmin         =>     ", inputs[:xmin])
        if(!haskey(inputs, :xmax))
            @warn " ::xmax is missing in .../IO/user_inputs.jl
           The default values :xmax=1.0 will be used."
            inputs[:xmax] = Float64(1.0)
        end
        println( " # :xmax         =>     ", inputs[:xmax])
        
        
        if (inputs[:nsd] > 1)
            if (!haskey(inputs, :npy))       
                @warn " :npy is missing in .../IO/user_inputs.
                The default value npy=2 will be used."
                inputs[:npy] = Int8(2)
                println( " # :npy          =>     ", inputs[:npy])
            else
                println( " # :npy          =>     ", inputs[:npy])
            end
            
            if(!haskey(inputs, :ymin))
                @warn " :ymin is missing in .../IO/user_inputs.jl
                The default value :ymin=-1.0 will be used."
                inputs[:ymin] = Float64(-1.0)
            end
            println( " # :ymin         =>     ", inputs[:ymin])
            if(!haskey(inputs, :ymin))
                @warn " ::ymax is missing in .../IO/user_inputs.
                The default values :ymax=1.0 will be used."
                inputs[:ymax] = Float64(1.0)
            end
            println( " # :ymax         =>     ", inputs[:ymax])

            
            if (inputs[:nsd] == 3)
                if (!haskey(inputs, :npz))
                    @warn " :npz is missing in .../IO/user_inputs.jl
                    The default value npz=2 will be used."
                    inputs[:npz] = Int8(2)
                    println( " # :npz          =>     ", inputs[:npz])
                else
                    println( " # :npz          =>     ", inputs[:npz])
                end
                
                if(!haskey(inputs, :zmin))
                    @warn " :zmin is missing in .../IO/user_inputs.jl
                    The default value :zmin=-1.0 will be used."
                    inputs[:zmin] = Float64(-1.0)
                end
                println( " # :zmin         =>     ", inputs[:zmin])
                if(!haskey(inputs, :zmax))
                    @warn " ::zmax is missing in .../IO/user_inputs.
                    The default values :zmax=1.0 will be used."
                    inputs[:zmax] = Float64(1.0)
                end                
                println( " # :zmax         =>     ", inputs[:zmax])
            end
        end
    end
    
    if error_flag == 1
        error(" JEXPRESSO inputs error: some necessary user inputs are not defined in IO/user_inputs.jl\n The program will stop now!")
    end

print(GREEN_FG(" # User inputs: ........................ DONE\n"))

#
# Correct quantities based on a hierarchy of input variables
#
if (inputs[:nsd] == 1)
    inputs[:npy] = 1
    inputs[:npz] = 1
elseif(inputs[:nsd] == 2)
    inputs[:npz] = 1
end

nvars::Int8 = 1
if (lowercase(inputs[:equation_set]) == "burgers")
    if(inputs[:nsd] == 1)
        nvars = 1
    elseif (inputs[:nsd] == 2)
        nvars = 2
    end
    println( " # nvars     ", nvars)
elseif (lowercase(inputs[:equation_set]) == "ns")
    if (inputs[:nsd] == 1)
        nvars = 3
    elseif(inputs[:nsd] == 2)
        nvars = 4
    elseif(inputs[:nsd] == 3)
        nvars == 5
    end
    println( " # nvars     ", nvars)
else
    println( "\n !!!! ERROR in user_inputs.jl: equation_set ", inputs[:equation_set], " is not coded!!!")
    println( "! Chose among: \n !!!!  [1] burgers1 \n")
    error()
end

return inputs, nvars
end

function mod_inputs_print_welcome()

    print(BLUE_FG(" #--------------------------------------------------------------------------------\n"))
    print(BLUE_FG(" # Welcome to ", RED_FG("jexpresso\n")))
    print(BLUE_FG(" # A Julia code to solve turbulence problems in the atmosphere\n"))
    print(BLUE_FG(" #--------------------------------------------------------------------------------\n"))

end

#inputs        = Dict{}()
#inputs, nvars = mod_inputs_user_inputs()

