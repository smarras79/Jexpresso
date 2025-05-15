# Tutorial: define the output variables:

This tutorial explains how to setup the array of output variables and extract them from the solution variables.

## Step 1: Go to your problem directory:

```bash
cd problems/equations/CompEuler/YOUR_PROBLEM
```

## Step 2: Open initialize.jl:

If not already done so, define the optional array of output variables `qoutvars` right after the solution array `q`.

The length of `qoutvars` and `qvars` do not have to be the same. In fact, `qoutvars` could possibly be way longer than `qvars`.

If `qoutvars` is not defined, the the default output variables coincide with `qvars`. Inside `initialize.jl` you may have:

```
    qvars    = ["ρ", "ρu", "ρv", "ρθ"]
    qoutvars = ["ρ", "u", "w", "θ", "p"]
```
 If `qoutvars` had not been previously used, remember to add it to the parameter argument of the call to `define_q` as shown below
    
    ```
    q = define_q(SD, 
                 mesh.nelem, 
                 mesh.npoin, 
                 mesh.ngl, 
                 qvars, 
                 TFloat, 
                 inputs[:backend]; 
                 neqs=length(qvars), 
                 qoutvars=qoutvars)
    ```

## Step 3: Extract the outout vartiables from the solutions:
Open `user_primitives,jl` and use `function user_out!` to derive the outpout vars from the unknown vector. 

Here an example to obtain some primitive quantities from the solution array `qvars=[ρ, ρu, ρv, ρθ]` if `qoutvars = ["ρ", "u", "w", "θ", "p"]` was defined in `initialize.jl` as explained above:

```
function user_uout!(uout, u, qe, )

    PhysConst = PhysicalConst{Float64}()

    uout[1] = u[1]       #ρ
    uout[2] = u[2]/u[1]  #u
    uout[3] = u[3]/u[1]  #v
    uout[4] = u[4]/u[1]  #θ

    p = perfectGasLaw_ρθtoP(PhysConst, ρ=uout[1], θ=uout[4])
    uout[end] = p
end

```

IMPORATANT NOTICE: the dimension of `uout` is that of `outvars` defined in `initialize.jl`. `uout` in `user_uout!` cannot be bigger than `qoutvars`!

Link to an example of [user_primitives.jl](https://github.com/smarras79/Jexpresso/blob/master/problems/equations/CompEuler/theta/user_primitives.jl).


