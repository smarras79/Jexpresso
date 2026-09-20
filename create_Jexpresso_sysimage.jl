using Pkg
using PackageCompiler
include("precompile_jexpresso.jl")
project_deps = Set(keys(Pkg.project().dependencies))
loaded = Set(k.name for (k,v) in Base.loaded_modules)
packages = sort(collect(intersect(project_deps, loaded)))
open("sysimage_packages.txt", "w") do f
           for p in packages; println(f, p); end
end
println("$(length(packages)) packages written")
packages = Symbol.(readlines("sysimage_packages.txt"))
create_sysimage(packages;
                sysimage_path = "jexpresso.so",
                precompile_execution_file = "precompile_jexpresso.jl")

