# Build the documentation (from the repository root):
#   julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#   julia --project=docs docs/make.jl
# The site is written to docs/build. Deployment (deploydocs) is not set up.
using Documenter
using ExponentialNoiseSDE

makedocs(
    sitename = "ExponentialNoiseSDE.jl",
    modules = [ExponentialNoiseSDE],
    pages = ["Home" => "index.md", "API" => "api.md"],
    checkdocs = :exports,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", edit_link = "main"),
)
