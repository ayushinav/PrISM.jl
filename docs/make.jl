using ProEM
using Documenter, DocumenterVitepress

include("pages.jl")
DocMeta.setdocmeta!(ProEM, :DocTestSetup, :(using ProEM); recursive=true)

makedocs(; authors="Abhinav Pratap Singh et al.",
    # sitename="ProEM.jl",
    # doctest=true,
    # clean=true,
    repo="github.com/ayushinav/ProEM.jl",
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="github.com/ayushinav/ProEM.jl", devurl="dev"),
    pages=pages,
    # draft=false,
    # source="src",
    # build="build"
    )

# DocumenterVitepress.deploydocs(; repo="github.com/ayushinav/ProEM.jl",
#     target=joinpath(@__DIR__, "build"), devbranch="main", push_preview=true)
