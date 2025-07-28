module Homework2
export de_casteljau, calculate_bezier_loop_area, calculate_bezier_loop_area_auto_detect, find_bezier_self_intersection

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

Izračuna določeni integral funkcije `f` na intervalu `[a, b]` z uporabo
Gauss-Legendre kvadrature z 20 vozlišči.

# Argumenti
- `f`: Funkcija, ki jo je treba integrirati. Pričakuje en argument `t::Float64`.
- `a`: Spodnja meja integracije.
- `b`: Zgornja meja integracije.

# Rezultat
- Vrednost določenega integrala.
"""
function integrate_gl_20(f::Function, a::Float64, b::Float64)::Float64
  total_integral = 0.0
  scale_factor = (b - a) / 2.0 # Skalacijski faktor za preslikavo intervala [-1, 1] na [a, b]

  for i in eachindex(GAUSS_LEGENDRE_NODES_20)
    node = GAUSS_LEGENDRE_NODES_20[i]
    weight = GAUSS_LEGENDRE_WEIGHTS_20[i]

    # Preslikava vozlišča (ξ) iz standardnega intervala [-1, 1] na interval [a, b]
    # t = ((b - a) / 2) * ξ + ((a + b) / 2)
    transformed_t = scale_factor * node + (a + b) / 2.0

    # Dodaj prispevek k integralu: utež * f(preslikano_vozlišče)
    total_integral += weight * f(transformed_t)
  end

  # Končni integral je skaliran s faktorjem (b - a) / 2
  return scale_factor * total_integral
end


"""
    bezier_loop_integrand(t::Float64, px::Vector{Float64}, py::Vector{Float64}, d_px::Vector{Float64}, d_py::Vector{Float64}) -> Float64

Izračuna vrednost integranda za izračun površine zanke Bézierove krivulje pri parametru `t`.
Integrand je \\(x(t)y'(t) - x'(t)y(t)\\).

# Argumenti
- `t`: Parameter na krivulji (običajno med 0.0 in 1.0).
- `px`: X-koordinate kontrolnih točk originalne Bézierove krivulje.
- `py`: Y-koordinate kontrolnih točk originalne Bézierove krivulje.
- `d_px`: X-koordinate kontrolnih točk derivata Bézierove krivulje (\\(x'(t)\\)).
- `d_py`: Y-koordinate kontrolnih točk derivata Bézierove krivulje (\\(y'(t)\\)).

# Rezultat
- Vrednost integranda \\(x(t)y'(t) - x'(t)y(t)\\) pri danem `t`.
"""
function bezier_loop_integrand(t::Float64, px::Vector{Float64}, py::Vector{Float64}, d_px::Vector{Float64}, d_py::Vector{Float64})::Float64
  # Izračunaj koordinate točke na krivulji pri t
  x_t = de_casteljau(px, t)
  y_t = de_casteljau(py, t)

  # Izračunaj koordinate odvoda krivulje pri t
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

  # Preverimo veljavnost vhodnih parametrov
  if !(0.0 <= t_start <= 1.0) || !(0.0 <= t_end <= 1.0) || t_start > t_end
    @warn "Neveljavni parametri t_start ali t_end. Integracija bo morda netočna."
    # Lahko bi vrnili 0.0 ali sprožili napako, odvisno od želenega obnašanja
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
    find_bezier_self_intersection(control_points::Vector{NTuple{2,Float64}}; num_samples_coarse::Int=1000, num_refinement_iterations::Int=100, refinement_samples::Int=10, distance_threshold_sq::Float64=1e-12, logging::Bool=true)

Finds the self-intersection parameters (t_start, t_end) of a Bézier curve using
a two-phase approach: a coarse search followed by iterative refinement.

# Arguments
- `control_points`: A vector of (x, y) tuples representing the control points.
- `num_samples_coarse`: Number of samples for the initial rough search.
- `num_refinement_iterations`: Number of iterations to refine the intersection point.
- `refinement_samples`: Number of samples to use within the shrinking window at each iteration.
- `distance_threshold_sq`: The squared distance threshold to accept a point as an intersection.
- `logging`: A boolean to enable or disable logging.

# Returns
A tuple `(t_start, t_end)` representing the parameter values at the point of self-intersection.
Throws an `ArgumentError` if no loop is found or if the refinement fails to converge.
"""
function find_bezier_self_intersection(
  control_points::Vector{NTuple{2,Float64}};
  num_samples_coarse::Int=1000,
  num_refinement_iterations::Int=100,
  refinement_samples::Int=10,
  distance_threshold_sq::Float64=1e-12,
  logging::Bool=true
)
  # ========== PHASE 1: COARSE SEARCH ==========
  t_values = range(0.0, stop=1.0, length=num_samples_coarse)
  points = [evaluate_bezier_2D(control_points, t) for t in t_values]
  min_index_separation = max(10, num_samples_coarse ÷ 20)

  min_dist_sq_coarse = Inf
  t_start_coarse = 0.0
  t_end_coarse = 0.0

  for i in 1:num_samples_coarse
    for j in (i+min_index_separation):num_samples_coarse
      p1 = points[i]
      p2 = points[j]
      dist_sq = (p1[1] - p2[1])^2 + (p1[2] - p2[2])^2
      if dist_sq < min_dist_sq_coarse
        min_dist_sq_coarse = dist_sq
        t_start_coarse = t_values[i]
        t_end_coarse = t_values[j]
      end
    end
  end

  if min_dist_sq_coarse > 1e-3 # Lenient threshold for coarse search
    throw(ArgumentError("No plausible loop found in coarse search. Min dist sq: $min_dist_sq_coarse"))
  end

  # ========== PHASE 2: ITERATIVE REFINEMENT ==========
  current_t_start = t_start_coarse
  current_t_end = t_end_coarse
  search_window = 2.0 / num_samples_coarse

  for iter in 1:num_refinement_iterations
    min_dist_sq_refinement = Inf
    best_t_start_in_iter = current_t_start
    best_t_end_in_iter = current_t_end

    t_start_range = range(clamp(current_t_start - search_window / 2, 0.0, 1.0), stop=clamp(current_t_start + search_window / 2, 0.0, 1.0), length=refinement_samples)
    t_end_range = range(clamp(current_t_end - search_window / 2, 0.0, 1.0), stop=clamp(current_t_end + search_window / 2, 0.0, 1.0), length=refinement_samples)

    for ts in t_start_range, te in t_end_range
      if abs(ts - te) < 1e-9
        continue
      end

      p1 = evaluate_bezier_2D(control_points, ts)
      p2 = evaluate_bezier_2D(control_points, te)
      dist_sq = (p1[1] - p2[1])^2 + (p1[2] - p2[2])^2

      if dist_sq < min_dist_sq_refinement
        min_dist_sq_refinement = dist_sq
        best_t_start_in_iter = ts
        best_t_end_in_iter = te
      end
    end

    current_t_start = best_t_start_in_iter
    current_t_end = best_t_end_in_iter
    search_window /= (refinement_samples / 2)
  end

  # ========== PHASE 3: FINAL CHECK ==========
  if current_t_start > current_t_end
    current_t_start, current_t_end = current_t_end, current_t_start
  end

  final_point_start = evaluate_bezier_2D(control_points, current_t_start)
  final_point_end = evaluate_bezier_2D(control_points, current_t_end)
  final_dist_sq = (final_point_start[1] - final_point_end[1])^2 + (final_point_start[2] - final_point_end[2])^2

  if final_dist_sq > distance_threshold_sq
    throw(ArgumentError("Refinement failed to converge. Min dist sq: $final_dist_sq > threshold $distance_threshold_sq"))
  end

  if logging
    println("Refined search complete. Final squared distance: $final_dist_sq")
    point_approx = ((final_point_start[1] + final_point_end[1]) / 2, (final_point_start[2] + final_point_end[2]) / 2)
    println("Detected point of intersection: $point_approx")
    println("Value of t_1:$current_t_start. Value of t_2:$current_t_end")
  end

  return (current_t_start, current_t_end)
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
