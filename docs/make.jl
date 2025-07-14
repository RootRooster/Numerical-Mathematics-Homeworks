using Documenter
using Homeworks

makedocs(
    sitename = "Homeworks",
    format = Documenter.HTML(),
    modules = [Homeworks]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
