module Homework3
using LinearAlgebra
using Plots

export rk4_step, calculate_distances, three_body_eom, simulate_trajectory_earth_moon_system, plot_trajectory, run_and_plot_simulation

# Physical Constants (using BigFloat for precision)
# https://en.wikipedia.org/wiki/Gravitational_constant
const G = 6.67430e-11 # Gravitational constant in m^3 kg^-1 s^-2

# https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
const M_EARTH = 5.9722e24 #Mass of Earth in kg

# https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
const M_MOON = 7.346e22 # Mass of Moon in kg

# distance between earth and moon centers
# https://en.wikipedia.org/wiki/Lunar_distance
const D_EM = 3.844e8 # Average Earth-Moon distance in m

# Derived dimensionless mass parameter mu
const MU_UNIT::Float64 = M_MOON / (M_EARTH + M_MOON)
const MU_BAR_UNIT::Float64 = 1 - MU_UNIT

"""
    calculate_distances(x::Float64, y::Float64, z::Float64)
Calculates the distance to the Earth (R) and the Moon (r) from a given point (x, y, z)
in a normalized coordinate system.

The function uses the following equations:
- R = sqrt((x + μ)^2 + y^2 + z^2)
- r = sqrt((x - μ̄)^2 + y^2 + z^2)

# Arguments
- `x::Float64`: The x-coordinate of the point.
- `y::Float64`: The y-coordinate of the point.
- `z::Float64`: The z-coordinate of the point.

# Returns
- `Tuple{Float64, Float64}`: A tuple containing the calculated distances (R, r).
"""
function calculate_distances(x::Float64, y::Float64, z::Float64)
  R = sqrt((x + MU_UNIT)^2 + y^2 + z^2)
  r = sqrt((x - MU_BAR_UNIT)^2 + y^2 + z^2)
  return (R, r)
end

"""
    rk4_step(f::Function, u::Vector{Float64}, h::Float64)
Performs a single step of the classic 4th-order Runge-Kutta method for solving ODEs.


# Arguments
- `f(u)`: The function defining the ODE system (e.g., `three_body_eom`).
- `u`: The current state vector.
- `h`: The step size.

# Returns
- The state vector at the next time step, approximated by the RK4 method.
"""
function rk4_step(f::Function, u::Vector{Float64}, h::Float64)
  # Note: .* is element wise multiplication
  k1 = h .* f(u)
  k2 = h .* f(u .+ 0.5 .* k1)
  k3 = h .* f(u .+ 0.5 .* k2)
  k4 = h .* f(u .+ k3)
  return u .+ (k1 .+ 2k2 .+ 2k3 .+ k4) ./ 6.0
end

"""
    three_body_eom(u::Vector{Float64}, mu::Float64)
Computes the equations of motion for the circular restricted three-body problem
in a non-dimensional, rotating frame.

This function is intended for use with an ODE solver (e.g., `rk4_step`).

# Arguments
- `u`: The 6-element state vector of the third body (e.g., a spacecraft):
    - `u[1:3]`: Position vector (x, y, z)
    - `u[4:6]`: Velocity vector (x_dot, y_dot, z_dot)

# Returns
- A 6-element vector representing the time derivative of the state vector:
    - `[x_dot, y_dot, z_dot, x_dot_dot, y_dot_dot, z_dot_dot]`

# Throws
- `ArgumentError`: If the input state vector `u` does not contain exactly 6 elements.
"""
function three_body_eom(u::Vector{Float64})
  # Check if the state vector has the correct number of elements
  if length(u) != 6
    throw(ArgumentError("Input state vector `u` must contain 6 elements (3 for position, 3 for velocity)."))
  end
  # Extract the values from the inital vector
  x, y, z = u[1], u[2], u[3] # These represent the position
  x_dot, y_dot, z_dot = u[4], u[5], u[6] # These represent the speed of the body

  # Pre-calculate common factors
  R, r = calculate_distances(x, y, z)
  mu_bar_over_R3 = MU_BAR_UNIT / (R^3)
  mu_over_r3 = MU_UNIT / (r^3)

  # Calculate the acceleration based on previous equations
  x_dot_dot = x + 2 * y_dot - mu_bar_over_R3 * (x + MU_UNIT) - mu_over_r3 * (x - MU_BAR_UNIT)
  # x_dot_dot = x - 2 * y_dot - mu_bar_over_R3 * (x + MU_UNIT) - mu_over_r3 * (x - MU_BAR_UNIT)
  y_dot_dot = y - 2 * x_dot - mu_bar_over_R3 * y - mu_over_r3 * y
  z_dot_dot = -mu_bar_over_R3 * z - mu_over_r3 * z

  # Return the derivatives of position and acceleration
  return [x_dot, y_dot, z_dot, x_dot_dot, y_dot_dot, z_dot_dot]
end

"""
    simulate_trajectory_earth_moon_system(u0::Vector{Float64}, t_span::Tuple{Float64, Float64}, h::Float64)
Simulates the trajectory of a spacecraft in the circular restricted three-body problem (Earth-Moon system).

# Arguments
- `u0::Vector{Float64}`: The initial 6-element state vector [x, y, z, x_dot, y_dot, z_dot].
- `t_span::Tuple{Float64, Float64}`: A tuple with the start and end time for the simulation (e.g., (0.0, 15.0)).
- `h::Float64`: The time step size for the RK4 integrator.

# Returns
- `Vector{Vector{Float64}}`: A vector of state vectors, representing the trajectory at each time step.
"""
function simulate_trajectory_earth_moon_system(u0::Vector{Float64}, t_span::Tuple{Float64,Float64}, h::Float64)
  # --- Argument Value Validation ---
  if length(u0) != 6
    throw(ArgumentError("Input state vector `u0` must contain 6 elements (3 for position, 3 for velocity)."))
  end
  if length(t_span) != 2
    throw(ArgumentError("The `t_span` parameter should be a tuple of two elements"))
  end
  if t_span[1] >= t_span[2]
    throw(ArgumentError("Start time in `t_span` must be less than the end time."))
  end
  if h <= 0.0
    throw(ArgumentError("Step size `h` must be positive."))
  end
  # --- End Validation ---

  # Define the equations of motion function to be passed to the solver
  eom_function = u -> three_body_eom(u)

  # Calculate the number of steps required for the simulation
  t_start, t_end = t_span
  num_steps = ceil(Int, (t_end - t_start) / h)

  # Pre-allocate a vector to store the entire trajectory for efficiency
  # The size is num_steps + 1 to include the initial state u0
  trajectory = Vector{Vector{Float64}}(undef, num_steps + 1)

  # Set the initial conditions
  trajectory[1] = u0
  u_current = u0

  # Main simulation loop
  for i in 1:num_steps
    # Calculate the next state using one RK4 step
    u_next = rk4_step(eom_function, u_current, h)

    # Store the new state and update the current state for the next iteration
    trajectory[i+1] = u_next
    u_current = u_next
  end

  return trajectory
end

"""
    plot_trajectory(trajectory::Vector{Vector{Float64}}, traj_type::String)
Creates and saves a 3D plot of the satellite's trajectory, with the
line color corresponding to the satellite's speed.

# Arguments
- `trajectory::Vector{Vector{Float64}}`: The trajectory history from the simulation.
- `traj_type::String`: A name for the trajectory type, used for the filename.
"""
function plot_trajectory(trajectory::Vector{Vector{Float64}}, traj_type::String)
  # --- Data Extraction ---
  # Extract position and calculate speed for each point in the trajectory.
  x_hist = [state[1] for state in trajectory]
  y_hist = [state[2] for state in trajectory]
  z_hist = [state[3] for state in trajectory]

  # Calculate the magnitude of the velocity vector (speed) at each state.
  speeds = [sqrt(state[4]^2 + state[5]^2 + state[6]^2) for state in trajectory]

  # --- Plotting Setup ---
  μ = MU_UNIT
  earth_pos = (-μ, 0, 0)
  moon_pos = (1 - μ, 0, 0)

  total_points = length(x_hist)
  sampling_factor = max(1, div(total_points, 2500))
  indices = 1:sampling_factor:total_points

  # --- Create the Plot ---
  p = plot3d(
    x_hist[indices], y_hist[indices], z_hist[indices],
    label="Trajectory",
    lc=:viridis,
    line_z=speeds[indices], # Color line by speed
    lw=1.5,
    xlabel="x (dimensionless)",
    ylabel="y (dimensionless)",
    zlabel="z (dimensionless)",
    title="Satellite Trajectory: $traj_type",
    legend=:topleft,
    camera=(30, 45),
    aspect_ratio=:equal,
    colorbar_title="Speed (dimensionless)" # label for the color bar
  )

  # Add Earth and Moon markers
  scatter!(p, [earth_pos[1]], [earth_pos[2]], [earth_pos[3]],
    markersize=8, markercolor=:dodgerblue, label="Earth")
  scatter!(p, [moon_pos[1]], [moon_pos[2]], [moon_pos[3]],
    markersize=4, markercolor=:grey, label="Moon")

  # --- Save the Figure ---
  filename = "trajectory_speed_$(traj_type).png"
  savefig(p, filename)
  println("Plot saved successfully to '$filename'")
end

"""
    run_and_plot_simulation(u0, t_span, h, traj_type)
Runs a full simulation and generates a plot of the resulting trajectory.

# Arguments
- `u0::Vector{Float64}`: The initial 6-element state vector.
- `t_span::Tuple{Float64, Float64}`: The simulation time span.
- `h::Float64`: The time step size for the RK4 integrator.
- `traj_type::String`: A name for the simulation, used for saving the plot.
"""
function run_and_plot_simulation(u0::Vector{Float64}, t_span::Tuple{Float64,Float64}, h::Float64, traj_type::String)
  println("Starting simulation for '$traj_type'...")
  trajectory = simulate_trajectory_earth_moon_system(u0, t_span, h)
  println("Simulation complete. Generating plot...")
  plot_trajectory(trajectory, traj_type)
end

end
