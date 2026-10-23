using Documenter
using PDBToolsBase
makedocs(
    modules = [PDBToolsBase],
    sitename = "PDBToolsBase.jl",
    repo = Remotes.GitHub("m3g", "PDBToolsBase.jl"),
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    # The examples in the docstrings use PDBTools.jl, and are tested in its documentation
    doctest = false,
    # Every public name must be documented in the API page
    checkdocs = :public,
    pages = [
        "Home" => "index.md",
        "Public API" => "api.md",
    ],
)
deploydocs(
    repo = "github.com/m3g/PDBToolsBase.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"],
)
