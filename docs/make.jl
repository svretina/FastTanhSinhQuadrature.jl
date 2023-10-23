using FastTanhSinhQuadrature
using Documenter

DocMeta.setdocmeta!(FastTanhSinhQuadrature, :DocTestSetup, :(using FastTanhSinhQuadrature); recursive=true)

makedocs(;
    modules=[FastTanhSinhQuadrature],
    authors="Stamatis Vretinaris",
    repo="https://github.com/svretina/FastTanhSinhQuadrature.jl/blob/{commit}{path}#{line}",
    sitename="FastTanhSinhQuadrature.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://svretina.github.io/FastTanhSinhQuadrature.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/svretina/FastTanhSinhQuadrature.jl",
    devbranch="master",
)
