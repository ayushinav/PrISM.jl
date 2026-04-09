using PrISM
using Documenter, DocumenterVitepress

include("pages.jl")
DocMeta.setdocmeta!(PrISM, :DocTestSetup, :(using PrISM); recursive=true)

makedocs(; authors="Abhinav Pratap Singh et al.",
    sitename="PrISM.jl",
    doctest=false,
    clean=true,
    repo="github.com/ayushinav/PrISM.jl",
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="github.com/ayushinav/PrISM.jl", devurl="dev", devbranch="main"),
    pages=pages,
    draft=false,
    source="src",
    build="build")

DocumenterVitepress.deploydocs(; repo="github.com/ayushinav/PrISM.jl",
    target=joinpath(@__DIR__, "build"), devbranch="main", push_preview=true)
