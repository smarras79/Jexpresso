using Documenter, Jexpresso

pages = [
  "Home" => "index.md",
 ]

makedocs(;
    modules=[Jexpresso],
    format=Documenter.HTML(),
    pages=pages,
    repo="https://github.com/smarras79/Jexpresso.jl/blob/{commit}{path}#L{line}",
    sitename="Jexpresso.jl",
    authors="S. Marras <smarras@njit.edu>, Y. Tissaoui <yt277@njit.edu>",
    # warnonly=true, # for debugging
)

deploydocs(;
    repo="github.com/smarras79/Jexpresso.jl",
)
