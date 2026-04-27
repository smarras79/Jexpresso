function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # SWE on the sphere -- stub case.
        # Equations, mesh, BCs, source, and initialization are not implemented yet.
        # This file exists so the case is registered and the problem boots.
        #---------------------------------------------------------------------------
        :case          => "swe_sphere",
        :SOL_VARS_TYPE => TOTAL(),
        :nimplemented  => true,
    )

    return inputs
end
