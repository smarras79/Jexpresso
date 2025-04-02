# This file is a part of project JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/NodeNumbering.jl/blob/master/LICENSE

using Documenter
using NodeNumbering

makedocs(
    modules = [NodeNumbering],
    sitename = "NodeNumbering.jl",
    format = :html,
    pages = [
             "Introduction" => "index.md",
             "Theory" => "theory.md",
             "API" => "api.md"
            ])
