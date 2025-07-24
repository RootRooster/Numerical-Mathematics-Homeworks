using Test
using Homework2
@testset "homework2" begin
  @testset "de_casteljau tests" begin

    function quadratic_bezier_test(points_quadratic::Vector{Float64}, t::Float64)
      return (1 - t)^2 * points_quadratic[1] + 2 * (1 - t) * t * points_quadratic[2] + t^2 * points_quadratic[3]
    end

    # Test case 1: Single point (degree 0)
    # The result should always be the point itself, regardless of t.
    @test de_casteljau([10.0], 0.5) == 10.0

    # Test case 2: Linear curve (degree 1)
    points_linear = [0.0, 10.0]
    @test de_casteljau(points_linear, 0.0) == 0.0
    @test de_casteljau(points_linear, 1.0) == 10.0
    @test de_casteljau(points_linear, 0.5) == 5.0 # Midpoint
    @test de_casteljau(points_linear, 0.2) ≈ 2.0 # (1-0.2)*0 + 0.2*10

    # Test case 3: Quadratic curve (degree 2)
    points_quadratic = [0.0, 10.0, 5.0]
    @test de_casteljau(points_quadratic, 0.0) == 0.0
    @test de_casteljau(points_quadratic, 1.0) == 5.0

    # Test against the explicit formula for a quadratic Bezier curve
    # B(t) = (1-t)^2*P0 + 2(1-t)t*P1 + t^2*P2
    t = 0.5
    expected_quadratic = quadratic_bezier_test(points_quadratic, t)
    @test de_casteljau(points_quadratic, t) ≈ expected_quadratic # Expected: 0.25*0 + 0.5*10 + 0.25*5 = 6.25

    t = 0.25
    expected_quadratic = quadratic_bezier_test(points_quadratic, t)
    @test de_casteljau(points_quadratic, t) ≈ expected_quadratic

    # Test case 4: Cubic curve (degree 3)
    points_cubic = [0.0, 10.0, 5.0, 20.0]
    @test de_casteljau(points_cubic, 0.0) == 0.0
    @test de_casteljau(points_cubic, 1.0) == 20.0

    # Test an intermediate point against the explicit formula
    t = 0.75
    p0, p1, p2, p3 = points_cubic
    expected_cubic = (1 - t)^3 * p0 + 3 * (1 - t)^2 * t * p1 + 3 * (1 - t) * t^2 * p2 + t^3 * p3
    @test de_casteljau(points_cubic, t) ≈ expected_cubic
  end

  @testset "calculate_bezier_loop_area tests (11-decimal precision)" begin
    # Test 1: Edge cases with no area
    #@test calculate_bezier_loop_area([]) == 0.0
    @test calculate_bezier_loop_area([(1.0, 1.0)]) == 0.0
    @test calculate_bezier_loop_area([(1.0, 1.0), (2.0, 2.0)]) == 0.0

    # Test 2: Collinear points that form a line segment (zero area)
    points_collinear = [(1.0, 1.0), (2.0, 2.0), (3.0, 3.0), (4.0, 4.0)]
    @test calculate_bezier_loop_area(points_collinear) ≈ 0.0 atol = 1e-11

    # Test 3: Quadratic curve (parabolic arc) with a known analytical area of -1/3.
    points_parabola = [(0.0, 0.0), (0.5, 1.0), (1.0, 0.0)]
    @test calculate_bezier_loop_area(points_parabola) ≈ -1 / 3 atol = 1e-11

    # Test 4: Reversing the parabola points should negate the signed area to +1/3.
    points_parabola_rev = reverse(points_parabola)
    @test calculate_bezier_loop_area(points_parabola_rev) ≈ 1 / 3 atol = 1e-11

    points_in_example = [
      (0.0, 0.0),
      (1.0, 1.0),
      (2.0, 3.0),
      (1.0, 4.0),
      (0.0, 4.0),
      (-1.0, 3.0),
      (0.0, 1.0),
      (1.0, 0.0)
    ]

    result = calculate_bezier_loop_area(points_in_example)
    @test result ≈ 1.966200466200 atol = 1e-11
  end
end
