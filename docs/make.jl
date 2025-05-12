#push!(LOAD_PATH,"../src/")

using Documenter, PrettyTables, Latexify
using Jexpresso

cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

pages = [
        "Home" => "index.md",
         "Performance" => "features/performance.md",
         "Best practices" => "features/best-practices.md",
         "User inputs" => "tutorials/user_inputs.md",
         "Tutorials" => Any["tutorials/running-jexpresso.md", 
                            "tutorials/newcase.md",
                            ],
        ]

makedocs(
    sitename = "Jexpresso.jl",
    authors = "Simone Marras, Yassine Tissaoui, Hang Wang",
    format = Documenter.HTML(
      size_threshold=nothing
    ),
    modules = [Jexpresso],
    pages = pages,
    doctest = false,
    warnonly = [:cross_references,:missing_docs],
    checkdocs = :exports,
    )
      
    deploydocs(
        repo = "github.com/smarras79/Jexpresso.jl.git",
    )