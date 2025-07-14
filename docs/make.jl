using Documenter
using Homeworks
using Homework1

makedocs(
  sitename="Homeworks Documentation",
  format=Documenter.HTML(
    edit_link="main",
    prettyurls=get(ENV, "CI", nothing) == "true",
  ),
  pages=[  # Your doc pages
    "Home" => "index.md",
    # Add more pages like "Homework 1" => "homework1.md"
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
  target="build",  # Where built HTML goes (default)
  push_preview=true  # Optional: Support PR previews
)

