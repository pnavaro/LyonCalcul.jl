using Documenter
using LyonCalcul

makedocs(modules=[LyonCalcul],
         doctest = false,
         format = Documenter.HTML(),
         sitename = "LyonCalcul.jl",
         pages = ["Documentation"    => "index.md",
                  "Talk" => "slides.md" ])

deploydocs(
    repo   = "github.com/pnavaro/LyonCalcul.jl.git",
 )

