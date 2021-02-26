using PlateReaderCore
using Documenter

DocMeta.setdocmeta!(PlateReaderCore, :DocTestSetup, :(using PlateReaderCore); recursive=true)

makedocs(;
    modules=[PlateReaderCore],
    authors="TAPO (Thomas Agersten Poulsen) <tapo@novozymes.com> and contributors",
    repo="https://github.com/tp2750/PlateReaderCore.jl/blob/{commit}{path}#{line}",
    sitename="PlateReaderCore.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tp2750.github.io/PlateReaderCore.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tp2750/PlateReaderCore.jl",
)
