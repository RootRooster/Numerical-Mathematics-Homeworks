using Test
using Homework3
using LinearAlgebra

# Set precision for tests to match module
setprecision(BigFloat, 128)

@testset "Homework3.jl Tests" begin
  @testset "RK4 Step Function Tests" begin
    @testset "1D Exponential Decay" begin
      u0 = [10.0]
      mu = 0.5
      h = 0.1
      u1_numeric = rk4_step(u -> -mu .* u, u0, h)
      t = h
      u1_exact = u0 .* exp(-mu * t)
      @test u1_numeric ≈ u1_exact atol = 1e-7
    end

    @testset "2D Simple Harmonic Oscillator" begin
      omega = 2.0
      mu = omega^2
      h = 0.01
      u0 = [1.0, 0.0]
      u1_numeric = rk4_step(u -> [u[2], -mu * u[1]], u0, h)
      t = h
      u1_exact = [cos(omega * t), -omega * sin(omega * t)]
      @test u1_numeric ≈ u1_exact atol = 1e-7
    end
  end

  @testset "Homework3: Distance Calculations" begin

    @testset "Origin Point" begin
      R, r = calculate_distances(0.0, 0.0, 0.0)
      @test R == Homework3.MU_UNIT
      @test r == Homework3.MU_BAR_UNIT
    end

    @testset "Earth's Location" begin
      R, r = calculate_distances(-Homework3.MU_UNIT, 0.0, 0.0)
      @test R ≈ 0.0 atol = 1e-15
      @test r ≈ 1.0
    end

    @testset "Moon's Location" begin
      R, r = calculate_distances(Homework3.MU_BAR_UNIT, 0.0, 0.0)
      @test R ≈ 1.0
      @test r ≈ 0.0 atol = 1e-15
    end
  end

  @testset "three_body_eom Tests" begin

    @testset "ArgumentError for incorrect vector size" begin
      @test_throws ArgumentError three_body_eom([1.0, 2.0, 3.0])
      @test_throws ArgumentError three_body_eom(Float64[])
      @test_throws ArgumentError three_body_eom([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    end

    # TODO: Maybe add Lagrange points equlibrium tests
  end

  @testset "simulate_trajectory tests" begin

    # A valid set of inputs to use for testing invalid cases
    valid_u0 = [1.1, 0.0, 0.0, 0.0, -0.1, 0.0]
    valid_t_span = (0.0, 1.0)
    valid_h = 0.1

    @testset "Argument validation" begin
      # Test for incorrect size of the initial state vector `u0`
      @test_throws ArgumentError simulate_trajectory_earth_moon_system([1.0, 2.0, 3.0], valid_t_span, valid_h)
      @test_throws ArgumentError simulate_trajectory_earth_moon_system([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0], valid_t_span, valid_h)

      # Test for invalid time span `t_span`
      @test_throws ArgumentError simulate_trajectory_earth_moon_system(valid_u0, (1.0, 0.0), valid_h) # End time before start time
      @test_throws ArgumentError simulate_trajectory_earth_moon_system(valid_u0, (1.0, 1.0), valid_h) # Start and end times are equal

      # Test for non-positive step size `h`
      @test_throws ArgumentError simulate_trajectory_earth_moon_system(valid_u0, valid_t_span, 0.0)
      @test_throws ArgumentError simulate_trajectory_earth_moon_system(valid_u0, valid_t_span, -0.01)
    end

    @testset "Successful simulation run" begin
      # Define parameters for a short, valid simulation
      u0 = [0.9, 0.0, 0.0, 0.0, 0.5, 0.0]
      t_span = (0.0, 2.0)
      h = 0.5

      # Expected number of steps = (2.0 - 0.0) / 0.5 = 4 steps.
      # Total trajectory points = 4 steps + 1 initial point = 5.
      expected_length = 5

      # Run the simulation
      trajectory = simulate_trajectory_earth_moon_system(u0, t_span, h)

      # Check the output
      @test trajectory isa Vector{Vector{Float64}}
      @test length(trajectory) == expected_length
      @test trajectory[1] == u0 # The first state should be the initial state
      @test length(trajectory[end]) == 6 # The final state should be a 6-element vector
    end

  end
end
