push!(LOAD_PATH,"../src/")

using Documenter
#using Jexpresso

makedocs(
    sitename="Jexpresso.jl",
    modules=[Jexpresso],
    format=Documenter.HTML(),
    pages = Any[
        "Home" => "index.md",
    ],
 )

deploydocs(;
           repo="github.com/smarras79/Jexpresso.jl",
)
