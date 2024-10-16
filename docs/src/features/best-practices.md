# Best practices to prevent performance deterioration

This section wants to give some non-comprehensive best-practice suggestions to prevent performance deterioration when new code is added or old code is modified.

Clearly, more thorough information on Julia performance can be found online, but here we report
those that we learned during the implementation of Jexpresso.


## Vectors vs tuples: 
do not abuse tuples. For example, the `user_inputs.jl` files that conbtain the setups for each problem case activate artificial diffusion as follows:
```
        :lvisc            => true,
        :ivisc_equations  => [1, 2, 3, 4],
        :μ                => [0.0, 20.0, 20.0, 60.0],
```
where the value of the diffusivity coefficient (`μ`) is given for each equation (`ivisc_equations`).

Notice that these are stored as arrays (`[...]`). While the code still works if we used tuples instead (`(...)`), performance would drasticially deterioriate due to unnecessary allocation at runtime.
Notice that tuples are perfectly fine to be used, as long as they are used where really needed.