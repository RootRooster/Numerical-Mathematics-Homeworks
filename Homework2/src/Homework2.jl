module Homework2
export de_casteljau, calculate_bezier_loop_area

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
    calculate_bezier_loop_area(control_points::Vector{NTuple{2, Float64}}) -> Float64

Calculates the signed area of a loop for a Bézier curve using Green's theorem.

# Arguments
- `control_points`: A vector of 2D tuples `(x, y)` representing the control polygon.

# Returns
- The signed area of the curve's loop.
"""
function calculate_bezier_loop_area(control_points::Vector{NTuple{2,Float64}})
  num_points = length(control_points)
  if num_points < 2
    return 0.0
  end

  degree = num_points - 1

  # Separate control points into x and y vectors
  px = [p[1] for p in control_points]
  py = [p[2] for p in control_points]

  # Calculate control points for the derivative curves
  d_px = degree * (px[2:end] - px[1:end-1])
  d_py = degree * (py[2:end] - py[1:end-1])

  total_integral = 0.0

  # Apply the Gauss-Legendre quadrature formula
  for i in eachindex(GAUSS_LEGENDRE_NODES_20)
    node = GAUSS_LEGENDRE_NODES_20[i]
    weight = GAUSS_LEGENDRE_WEIGHTS_20[i]

    # The formula works on the interval [-1, 1], but our parameter `t` is in [0, 1].
    # We transform the node to our interval.
    t = 0.5 * (node + 1.0)

    # Evaluate curve and its derivative at t
    x_t = de_casteljau(px, t)
    y_t = de_casteljau(py, t)
    dx_t = de_casteljau(d_px, t)
    dy_t = de_casteljau(d_py, t)

    # Calculate the function we want to integrate
    integrand = x_t * dy_t - dx_t * y_t
    total_integral += weight * integrand
  end

  # The final area calculation requires two factors of 1/2:
  # 1. One from the area formula itself: P = (1/2) * integral(...)
  # 2. One from changing the integration interval from [0,1] to [-1,1] for the quadrature.
  # Total factor is (1/2) * (1/2) = 0.25.
  area = 0.25 * total_integral

  return area
end


end # module Homework2
