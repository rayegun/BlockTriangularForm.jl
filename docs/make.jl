using BlockTriangularForm
using Documenter

DocMeta.setdocmeta!(BlockTriangularForm, :DocTestSetup, :(using BlockTriangularForm); recursive=true)

makedocs(;
    modules=[BlockTriangularForm],
    authors="Wimmerer <kimmerer@mit.edu> and contributors",
    repo="https://github.com/"Wimmerer"/BlockTriangularForm.jl/blob/{commit}{path}#{line}",
    sitename="BlockTriangularForm.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://"Wimmerer".github.io/BlockTriangularForm.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/"Wimmerer"/BlockTriangularForm.jl",
    devbranch="main",
)
