using Test
using Pkg

const ROOT = abspath(joinpath(@__DIR__, ".."))

@testset "All Homeworks Tests" begin
  @testset "Main Homeworks Project" begin
    @test true
  end

  # --- Test Sub-Projects ---
  # Define the list of sub-project folders to test
  sub_projects = ["Homework1", "Homework2", "Homework3"]

  for project_name in sub_projects
    @testset "$project_name Tests" begin
      project_path = joinpath(ROOT, project_name)

      if !isdir(project_path)
        @error "Sub-project directory not found: $project_path"
        @test false # Mark this testset as failed
        continue
      end

      # Temporarily change to the sub-project's directory
      original_dir = pwd()
      try
        cd(project_path)

        # Activate the sub-project's specific environment
        Pkg.activate(".")

        # Install its dependencies
        Pkg.instantiate()

        # Run its tests by including its runtests.jl file
        include(joinpath(project_path, "test", "runtests.jl"))
      finally
        # Always return to the original directory and reactivate the main environment
        cd(original_dir)
        Pkg.activate(ROOT)
      end
    end
  end
end
