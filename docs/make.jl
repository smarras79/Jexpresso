push!(LOAD_PATH,"../src/")

using Documenter, Jexpresso

makedocs(
    sitename="Jexpresso.jl",
    modules=[Jexpresso],
    format=Documenter.HTML(),
    pages = Any[
        "Index" => "index1.md",
        "Jexpresso" => "Jexpresso.md",
        "`_build_rhs!`" => "Build rhs.md",
        "`_expansion_inviscid!`" => "Expansion inviscid.md",
        "Performance Comparison" => "Performance Comparison.md",
    ],
 )

 deploydocs(
    repo = "github.com/smarras79/Jexpresso.jl.git",

)
