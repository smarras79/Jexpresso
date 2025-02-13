#push!(LOAD_PATH,"../src/")

using Documenter, PrettyTables, Latexify
using Jexpresso

cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(sitename = "Jexpresso.jl",
         authors = "Simone Marras, Yassine Tissaoui, Hang Wang",
         format = Documenter.HTML(assets = ["assets/favicon.ico"],
                                  size_threshold=nothing
                                  ),
         
         modules = [Jexpresso],
         doctest = false,
         clean = true,
         pages = pages,
         warnonly = [:cross_references,:missing_docs],
         checkdocs = :exports,
         )

deploydocs(repo="github.com/smarras79/Jexpresso.jl.git";
           push_preview = true)
