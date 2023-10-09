push!(LOAD_PATH,"../src/")

using Documenter, Jexpresso

makedocs(
    sitename="Jexpresso.jl",
    modules=[Jexpresso],
    format=Documenter.HTML(),
    pages = Any[
        "Home" => "index.md",
    ],
    #repo="https://github.com/smarras79/Jexpresso.jl/blob/{commit}{path}#L{line}",
    #authors="S. Marras <smarras@njit.edu>, Y. Tissaoui <yt277@njit.edu>"
    # warnonly=true, # for debugging
)

deploydocs(;
    repo="https://github.com/smarras79/Jexpresso.jl/blob/{commit}{path}#L{line}",
#           repo="github.com/smarras79/Jexpresso.jl",
)
