#push!(LOAD_PATH,"../src/")

using Documenter, Jexpresso, PrettyTables, Latexify

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(
    sitename = "Jexpresso.jl",
    format = Documenter.HTML(assets = ["assets/favicon.ico"],
                             canonical = "https://smarras79.github.io/Jexpresso/dev/",
                             size_threshold=nothing
                             ),
    authors = "Simone Marras, Yassine Tissaoui, Hang Wang",
    modules = [Jexpresso],
    #linkcheck = true,
    doctest = false, clean = true,
    pages = pages,
    warnonly = [:cross_references,:missing_docs],
    checkdocs = :exports,
    )

deploydocs(repo="github.com/smarras79/Jexpresso.jl")
