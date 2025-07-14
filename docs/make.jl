using Documenter
using Homeworks
using Homework1

makedocs(
  sitename="Homeworks Documentation",
  pages=[
    "Overview" => "index.md",
  ],
  modules=[Homeworks, Homework1]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
deploydocs(
  repo="github.com/RootRooster/Numerical-Mathematics-Homeworks",  # Your repo URL (no https://)
  devbranch="main",  # Branch where your source code lives
)

