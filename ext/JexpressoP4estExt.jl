"""
Jexpresso's AMR code path: everything that calls into GridapP4est /
P4est_wrapper. Loaded automatically when GridapP4est is loaded alongside
Jexpresso; see src/kernel/Adaptivity/p4est_hooks.jl for the hooks it fills
and for how Jexpresso loads it.
"""
module JexpressoP4estExt

import Jexpresso as JX
import GridapP4est

# The implementation is shared with the script-mode loader
# (`julia src/Jexpresso.jl ...` is not a package, so no extension), which
# includes the same file into a module with the same `JX`/`GridapP4est`
# bindings.
include(joinpath("JexpressoP4estExt", "impl.jl"))

function __init__()
    JX._AMR_LOADED[] = true
end

end # module
