using Documenter
using LyonCalcul

makedocs(modules=[LyonCalcul],
         doctest = false,
         format = Documenter.HTML(),
         sitename = "LyonCalcul.jl",
         pages = ["Documentation"    => "index.md"])

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/pnavaro/LyonCalcul.jl.git",
 )

