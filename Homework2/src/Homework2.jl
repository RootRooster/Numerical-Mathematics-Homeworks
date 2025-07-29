module Homework2
using LinearAlgebra
export de_casteljau, calculate_bezier_loop_area, calculate_bezier_loop_area_auto_detect, find_bezier_self_intersection, bezier_loop_integrand, integrate_gl_20, evaluate_bezier_2D, calculate_bezier_loop_area, derivative_control_points

# Pre-computed 20-point Gauss-Legendre quadrature nodes and weights for the interval [-1, 1].
# Source: https://dlmf.nist.gov/3.5
# This provides extremely high precision, sufficient for integrating polynomials up to degree 39 exactly.
const GAUSS_LEGENDRE_NODES_20 = [
  -0.9931285991850949, -0.9639719272779138, -0.9122344282513259, -0.8391169718222188,
  -0.7463319064601508, -0.6360536807265150, -0.5108670019508271, -0.3737060887154195,
  -0.2277858511416451, -0.0765265211334973, 0.0765265211334973, 0.2277858511416451,
  0.3737060887154195, 0.5108670019508271, 0.6360536807265150, 0.7463319064601508,
  0.8391169718222188, 0.9122344282513259, 0.9639719272779138, 0.9931285991850949
]

const GAUSS_LEGENDRE_WEIGHTS_20 = [
  0.0176140071391521, 0.0406014298003869, 0.0626720483341091, 0.0832767415767048,
  0.1019301198172404, 0.1181945319615184, 0.1316886384491766, 0.1420961093183821,
  0.1491729864726037, 0.1527533871307259, 0.1527533871307259, 0.1491729864726037,
  0.1420961093183821, 0.1316886384491766, 0.1181945319615184, 0.1019301198172404,
  0.0832767415767048, 0.0626720483341091, 0.0406014298003869, 0.0176140071391521
]

"""
    de_casteljau(points::Vector{Float64}, t::Float64) -> Float64

Evaluates a 1D Bézier curve at parameter `t` using de Casteljau's algorithm.

# Arguments
- `points`: A vector of control point coordinates (e.g., all x-coordinates).
- `t`: The parameter value, typically in [0, 1].

# Returns
- The value of the curve at `t`.
"""
function de_casteljau(points::Vector{Float64}, t::Float64)
  n = length(points)
  if n == 1
    return points[1]
  end

  # Create a mutable copy of the points to work with
  coeffs = copy(points)

  # Iteratively apply linear interpolation
  for i in 1:(n-1)
    for j in 1:(n-i)
      coeffs[j] = (1.0 - t) * coeffs[j] + t * coeffs[j+1]
    end
  end

  return coeffs[1]
end

"""
    evaluate_bezier(control_points::Vector{NTuple{2,Float64}}, t::Float64)

Evaluates a 2D Bézier curve at parameter `t` using de Casteljau's algorithm.

# Arguments
- `control_points`: A vector of control point coordinates (e.g., all x-coordinates).
- `t`: The parameter value, typically in [0, 1].

# Returns
- The value (x,y) of the curve at `t`.
"""
function evaluate_bezier_2D(control_points::Vector{NTuple{2,Float64}}, t::Float64)
  n = length(control_points) - 1
  if n < 0
    return (0.0, 0.0)
  end
  # Extract the x and y coordinates into separate vectors
  x_points = [p[1] for p in control_points]
  y_points = [p[2] for p in control_points]

  # Calculate the interpolated x and y coordinates using de_casteljau
  x = de_casteljau(x_points, t)
  y = de_casteljau(y_points, t)

  return (x, y)
end

"""
    integrate_gl_20(f::Function, a::Float64, b::Float64) -> Float64

Calculates the definite integral of a function `f` over the interval `[a, b]`
using 20-point Gauss-Legendre quadrature.

# Arguments
- `f`: The function to integrate. It should accept a single `Float64` argument.
- `a`: The lower bound of integration.
- `b`: The upper bound of integration.

# Returns
- The approximate value of the definite integral.
"""
function integrate_gl_20(f::Function, a::Float64, b::Float64)::Float64
  total_integral = 0.0
  # Scaling factor to map the standard interval [-1, 1] to [a, b]
  scale_factor = (b - a) / 2.0

  for i in eachindex(GAUSS_LEGENDRE_NODES_20)
    node = GAUSS_LEGENDRE_NODES_20[i]
    weight = GAUSS_LEGENDRE_WEIGHTS_20[i]

    # Map the node (ξ) from the standard interval [-1, 1] to the target interval [a, b]
    # The transformation is: t = ((b - a) / 2) * ξ + ((a + b) / 2)
    transformed_t = scale_factor * node + (a + b) / 2.0

    # Add the contribution to the integral: weight * f(transformed_node)
    total_integral += weight * f(transformed_t)
  end

  # The final integral must be scaled by the factor (b - a) / 2
  return scale_factor * total_integral
end

"""
    bezier_loop_integrand(t::Float64, px::Vector{Float64}, py::Vector{Float64}, d_px::Vector{Float64}, d_py::Vector{Float64}) -> Float64

Calculates the value of the integrand used to compute the area of a Bézier curve's loop at a given parameter `t`.

The area of a closed loop can be computed using Green's Theorem, where the integrand is given by \\(x(t)y'(t) - x'(t)y(t)\\).

# Arguments
- `t`: The parameter on the curve (typically between 0.0 and 1.0).
- `px`: The x-coordinates of the control points for the original Bézier curve.
- `py`: The y-coordinates of the control points for the original Bézier curve.
- `d_px`: The x-coordinates of the control points for the derivative Bézier curve (representing \\(x'(t)\\)).
- `d_py`: The y-coordinates of the control points for the derivative Bézier curve (representing \\(y'(t)\\)).

# Returns
- The value of the integrand \\(x(t)y'(t) - x'(t)y(t)\\) at the given `t`.
"""
function bezier_loop_integrand(t::Float64, px::Vector{Float64}, py::Vector{Float64}, d_px::Vector{Float64}, d_py::Vector{Float64})::Float64
  # Calculate the value at t
  x_t = de_casteljau(px, t)
  y_t = de_casteljau(py, t)

  # Calculate the derivative value at t
  dx_t = de_casteljau(d_px, t)
  dy_t = de_casteljau(d_py, t)

  return x_t * dy_t - dx_t * y_t
end

"""
    calculate_bezier_loop_area(control_points::Vector{NTuple{2,Float64},t_start::Float64 = 0.0,t_end::Float64 = 1.0}) -> Float64

Calculates the signed area of a loop for a Bézier curve using Green's theorem.

# Arguments
- `control_points`: A vector of 2D tuples `(x, y)` representing the control polygon.
- `t_start`: Start parameter (t-value) for integration. Default is 0.0.
- `t_end`: Končni parameter (t-vrednost) za integracijo. Default is 1.0.

# Returns
- The signed area of the curve's loop.
"""
function calculate_bezier_loop_area(control_points::Vector{NTuple{2,Float64}}; t_start::Float64=0.0, t_end::Float64=1.0)
  num_points = length(control_points)
  if num_points < 2
    return 0.0
  end

  if !(0.0 <= t_start <= 1.0) || !(0.0 <= t_end <= 1.0) || t_start > t_end
    @warn "Neveljavni parametri t_start ali t_end. Integracija bo morda netočna."
  end

  degree = num_points - 1

  # Separate control points into x and y vectors
  px = [p[1] for p in control_points]
  py = [p[2] for p in control_points]

  # Calculate control points for the derivative curves
  d_px = degree * (px[2:end] - px[1:end-1])
  d_py = degree * (py[2:end] - py[1:end-1])

  # Define an anonymous function for the integrand that captures the calculated control points of the derivatives.
  # This function will be passed to the integration function.
  current_bezier_integrand = t -> bezier_loop_integrand(t, px, py, d_px, d_py)

  # Calculate the definite integral of the integrand from 0.0 to 1.0. 
  # According to Green's theorem, the area A = 0.5 * integral(x dy - y dx).
  # The function integrate_gl_20 returns the value of the integral, so we multiply it by 0.5.
  loop_integral_value = integrate_gl_20(current_bezier_integrand, t_start, t_end)
  return 0.5 * loop_integral_value
end

"""
    derivative_control_points(points::Vector{NTuple{2, Float64}})

Calculates the control points for the derivative of a Bézier curve.

The derivative of a degree (n) Bézier curve is a Bézier curve of degree (n-1). This function computes the control points for that derivative curve.

# Arguments
- `points`: A vector of 2D tuples representing the control points of the original curve.

# Returns
- A new vector of control points for the derivative curve.
"""
function derivative_control_points(points::Vector{NTuple{2,Float64}})
  n = length(points) - 1
  return [(n * (points[i+1][1] - points[i][1]), n * (points[i+1][2] - points[i][2])) for i in 1:n]
end

"""
    find_bezier_self_intersection(control_points::Vector{NTuple{2,Float64}}; kwargs...) -> Tuple{Float64, Float64}

Finds the self-intersection parameters (t_start, t_end) of a Bézier curve.
This uses a two-phase approach: a coarse search to find an initial guess,
then refines it with Newton's method for solving systems of nonlinear equations.

# Arguments
- `control_points`: A vector of `(x, y)` tuples representing the control points.
- `num_samples_coarse::Int=1000`: Number of samples for the initial coarse search.
- `max_newton_iter::Int=20`: Maximum number of iterations for Newton's method.
- `tolerance::Float64=1e-12`: The convergence tolerance for Newton's method (based on the norm of the step).
- `logging::Bool=true`: A boolean to enable or disable logging output.

# Returns
- A sorted tuple `(t_start, t_end)` representing the parameter values at the point of self-intersection.

# Throws
- `ArgumentError`: If no plausible self-intersection is found in the coarse search.
- `ArgumentError`: If Newton's method fails to converge or if the Jacobian matrix becomes singular.
"""
function find_bezier_self_intersection(
  control_points::Vector{NTuple{2,Float64}};
  num_samples_coarse::Int=1000,
  max_newton_iter::Int=20,
  tolerance::Float64=1e-12,
  logging::Bool=true
)
  # STEP 1: Coarse search for a good initial approximation
  samples = [evaluate_bezier_2D(control_points, t) for t in range(0.0, 1.0, length=num_samples_coarse)]
  t_values = collect(range(0.0, 1.0, length=num_samples_coarse))

  min_dist_sq = Inf
  t1_coarse, t2_coarse = -1.0, -1.0

  for i in 1:num_samples_coarse
    for j in (i+1):num_samples_coarse
      # Avoid comparing adjacent points on the curve
      if abs(t_values[i] - t_values[j]) < 0.1
        continue
      end

      dist_sq = (samples[i][1] - samples[j][1])^2 + (samples[i][2] - samples[j][2])^2
      if dist_sq < min_dist_sq
        min_dist_sq = dist_sq
        t1_coarse, t2_coarse = t_values[i], t_values[j]
      end
    end
  end

  if t1_coarse == -1.0
    throw(ArgumentError("No self-intersection found in coarse search."))
  end

  if logging
    println("Coarse search found initial guess: t1 ≈ $t1_coarse, t2 ≈ $t2_coarse")
  end

  # STEP 2: Refined search using Newton's method
  t = [t1_coarse, t2_coarse]
  d_control_points = derivative_control_points(control_points)

  for iter in 1:max_newton_iter
    # Define the function F(t) = B(t₁) - B(t₂) which we want to find the root of
    val_t1 = evaluate_bezier_2D(control_points, t[1])
    val_t2 = evaluate_bezier_2D(control_points, t[2])
    F = [val_t1[1] - val_t2[1], val_t1[2] - val_t2[2]]

    # Jacobian matrix J_F(t) = [B'(t₁) | -B'(t₂)]
    d_val_t1 = evaluate_bezier_2D(d_control_points, t[1])
    d_val_t2 = evaluate_bezier_2D(d_control_points, t[2])
    J = [d_val_t1[1] -d_val_t2[1];
      d_val_t1[2] -d_val_t2[2]]

    # Check for matrix singularity
    if abs(det(J)) < 1e-12
      throw(ArgumentError("Jacobian matrix is singular. Newton's method failed."))
    end

    # Solve the system J * delta_t = -F for the Newton step
    delta_t = J \ (-F)
    t += delta_t

    # Stopping condition
    if norm(delta_t) < tolerance
      if logging
        intersection_point = evaluate_bezier_2D(control_points, t[1])
        println("Newton's method converged in $iter iterations.")
        println("Detected intersection point: $intersection_point")
        println("Final parameters: t1 = $(t[1]), t2 = $(t[2])")
      end
      return sort(t) # Return the sorted pair of t-values
    end
  end

  throw(ErrorException("Newton's method did not converge within $max_newton_iter iterations."))
end

"""
    calculate_bezier_loop_area_auto_detect(control_points::Vector{NTuple{2,Float64}}; kwargs...)

Finds the self-intersection of a Bézier curve and calculates the area of the loop.
This function first calls `find_bezier_self_intersection` to determine the
intersection parameters `t_start` and `t_end`.

# Arguments
- `control_points`: A vector of (x, y) tuples representing the control points.
- `...`: Keyword arguments to be passed to `find_bezier_self_intersection`.

# Returns
The area of the Bézier curve loop.
"""
function calculate_bezier_loop_area_auto_detect(
  control_points::Vector{NTuple{2,Float64}};
  kwargs...
)
  # Find the self-intersection parameters t_start and t_end
  t_start, t_end = find_bezier_self_intersection(control_points; kwargs...)

  # Calculate the area of the loop defined by these parameters
  return calculate_bezier_loop_area(
    control_points,
    t_start=t_start,
    t_end=t_end
  )
end


end # module Homework2
