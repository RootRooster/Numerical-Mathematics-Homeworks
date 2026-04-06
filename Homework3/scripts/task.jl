using Homework3

println("Earth position: $(Homework3.MU_UNIT)")
println("Moon position: $(Homework3.MU_BAR_UNIT)")

# --------------- SIMPLE ORBIT SCENARIO ---------------  
println("\n--- Running a Simple Orbit---")
# Initial conditions for a simple flyby
u0_flyby = [0.2, 0.0, 0.0, 0.0, 1.938, 0.00]
t_span_flyby = (0.0, 1.0)
h_flyby = 0.005

# Run the second simulation
run_and_plot_simulation(u0_flyby, t_span_flyby, h_flyby, "simple_orbit")



# --------------- SIMPLE FREE RETURN TRAJECTORY ---------------  
println("\n--- Running a Free Return trajectory ---")
u0_flyby = [0.2, 0.0, 0.0, 1.954, 1.78, 0.0]

t_span_flyby = (0.0, 10.0)
h_flyby = 0.005

# Run the second simulation
run_and_plot_simulation(u0_flyby, t_span_flyby, h_flyby, "free_return_trajectory")


