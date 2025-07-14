using Documenter
using Homeworks
using Homework1

makedocs(
  sitename="Homeworks Documentation",
  format=Documenter.HTML(
    edit_link="homework_1"
  ),
  modules=[Homeworks, Homework1]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
