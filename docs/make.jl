using CBChPT
using Documenter

DocMeta.setdocmeta!(CBChPT, :DocTestSetup, :(using CBChPT); recursive=true)

makedocs(;
    modules=[CBChPT],
    authors="Zejian Zhuang zejian.zhuang@uv.es",
    sitename="CBChPT.jl",
    format=Documenter.HTML(;
        canonical="https://zejianzhuang-uv.github.io/CBChPT.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/zejianzhuang-uv/CBChPT.jl",
    devbranch="main",
)
