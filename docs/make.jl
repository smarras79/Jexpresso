push!(LOAD_PATH,"../src/")

using Documenter, Jexpresso

makedocs(
    sitename="Jexpresso.jl",
    modules=[Jexpresso],
    format=Documenter.HTML(),
    pages = Any[
        "Home" => "index.md",
        "Jexpresso" => "Jexpresso.md",
        "_build_rhs!" => "Build rhs.md",
        "_expansion_inviscid!" => "Expansion inviscid.md",
    ],
 )

 deploydocs(
    repo = "github.com/smarras79/Jexpresso.jl.git",
    target = "gh-pages",
)
