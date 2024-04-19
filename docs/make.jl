push!(LOAD_PATH,"../src/")

using Documenter
#using Jexpresso

pages = [
  "Home" => "index.md",
]

makedocs(
  sitename = "Jexpresso.jl",
  format = Documenter.HTML(
    size_threshold=nothing
  ),
  modules = [Jexpresso],
  pages = pages,
  doctest = false,
  warnonly = [:cross_references,:missing_docs],
  checkdocs = :exports,
)

deploydocs(;
           repo="github.com/smarras79/Jexpresso.jl",
)
