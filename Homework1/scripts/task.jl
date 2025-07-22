# License: © 2025 Martin Vuk. All rights reserved. Explicit permission by the author to use the code in this solution for the homework.
# I have obtained verbal consent from the author to utilize this code for the purpose of this project.

using Graphs
using Plots
using GraphRecipes
using Homework1

"""
    G = krožna_lestev(n)

Ustvari graf krožna lestev z `2n` točkami.
"""
function krožna_lestev(n)
  G = SimpleGraph(2 * n)
  # prvi cikel
  for i = 1:n-1
    add_edge!(G, i, i + 1)
  end
  add_edge!(G, 1, n)
  # drugi cikel
  for i = n+1:2n-1
    add_edge!(G, i, i + 1)
  end
  add_edge!(G, n + 1, 2n)
  # povezave med obema cikloma
  for i = 1:n
    add_edge!(G, i, i + n)
  end
  return G
end

"""
    A = matrika(G::AbstractGraph, sprem)

Poišči matriko sistema linearnih enačb za vložitev grafa `G` s fizikalno metodo.
Argument `sprem` je vektor vozlišč grafa, ki nimajo določenih koordinat.
Indeksi v matriki `A` ustrezajo vozliščem v istem vrstnem redu,
kot nastopajo v argumentu `sprem`.
Vrne matriko tipa RedkaMatrika
"""
function Rmatrikafy(G::AbstractGraph, sprem)
  # preslikava med vozlišči in indeksi v matriki
  v_to_i = Dict([sprem[i] => i for i in eachindex(sprem)])
  n = length(sprem)
  V_vecs = [Vector{Float64}() for _ in 1:n]
  I_vecs = [Vector{Int}() for _ in 1:n]
  for i = 1:n
    vertex = sprem[i]
    sosedi = neighbors(G, vertex)
    for vertex2 in sosedi
      if haskey(v_to_i, vertex2)
        j = v_to_i[vertex2]
        push!(V_vecs[i], 1.0)
        push!(I_vecs[i], j)
      end
    end
    push!(V_vecs[i], -length(sosedi))
    push!(I_vecs[i], i)
  end
  return RedkaMatrika(V_vecs, I_vecs)
end

"""
    b = desne_strani(G::AbstractGraph, sprem, koordinate)

Poišči desne strani sistema linearnih enačb za eno koordinato vložitve grafa `G`
s fizikalno metodo. Argument `sprem` je vektor vozlišč grafa, ki nimajo
določenih koordinat. Argument `koordinate` vsebuje eno koordinato za vsa
vozlišča grafa. Metoda uporabi le koordinato vozlišč, ki so pritrjena.
Indeksi v vektorju `b` ustrezajo vozliščem v istem vrstnem redu,
kot nastopajo v argumentu `sprem`.
"""
function desne_strani(G::AbstractGraph, sprem, koordinate)
  set = Set(sprem)
  m = length(sprem)
  b = zeros(m)
  for i = 1:m
    v = sprem[i]
    for v2 in neighbors(G, v)
      if !(v2 in set) # dodamo le točke, ki so fiksirane
        b[i] -= koordinate[v2]
      end
    end
  end
  return b
end

"""
  negative(A::RedkaMatrika)

Vrne matriko z obratno predznačenimi elementi.
"""
function negative(A::RedkaMatrika)
  I = deepcopy(A.I)
  V = deepcopy(A.V)
  n = length(A.V)
  for row_values in V
    row_values .= -row_values
  end
  return RedkaMatrika(V, I)
end

"""
    vloži!(G::AbstractGraph, fix, točke)

Poišči vložitev grafa `G` v prostor s fizikalno metodo. Argument `fix` vsebuje
vektor vozlišč grafa, ki imajo določene koordinate. Argument `točke` je
začetna vložitev grafa. Koordinate vozlišč, ki niso pritrjena, bodo nadomeščene
z novimi koordinatami.

Metoda ne vrne ničesar, ampak zapiše izračunane koordinate v matriko `točke`.
"""
function vloži!(G::AbstractGraph, fix, točke)
  sprem = setdiff(vertices(G), fix)
  dim, _ = size(točke)
  A = Rmatrikafy(G, sprem)
  neg_A = negative(A)
  for k = 1:dim
    b = desne_strani(G, sprem, točke[k, :])
    x0 = zeros(length(b))
    x, _ = sor(neg_A, -b, x0, 0.2) # matrika A je negativno definitna
    točke[k, sprem] = x
  end
end

#
# Uporaba metode za vložitev grafa v ravnino s fizikalno metodo.
#
println("Vlagam graf v ravnino.")
G = krožna_lestev(8)
graphplot(G, curves=false)
t = range(0, 2pi, 9)[1:end-1]
x = cos.(t)
y = sin.(t)
točke = hcat(hcat(x, y)', zeros(2, 8))
# funkcija hcat zloži vektorje po stolpcih v matriko
fix = 1:8
vloži!(G, fix, točke)

final_plot = graphplot(G, x=točke[1, :], y=točke[2, :], curves=false, aspect_ratio=:equal, title="Circular Ladder Graph (n=8)")
savefig(final_plot, "docs/src/assets/krozna_lestev_plot.png")
println("Plot shranjen v 'docs/src/assets/krozna_lestev_plot.png'")


#
# Iskanje optimalenga parametra omega
#

# Pripravimo sistem iz primera
sprem = setdiff(vertices(G), fix)
A = Rmatrikafy(G, sprem)
# Vzamemo prvo koordinato (x)
b = desne_strani(G, sprem, točke[1, :])
x0 = zeros(length(b))

# Testiraj različne omega vrednosti
omega_values = 0.1:0.1:1.9
iterations_list = []

neg_A = negative(A)
for omega in omega_values
  try
    local x, iterations = sor(neg_A, -b, x0, omega)
    push!(iterations_list, iterations)
  catch
    push!(iterations_list, NaN)  # če ne konvergira
  end
end

# Poišči optimalno omega
valid_indices = .!isnan.(iterations_list)
if any(valid_indices)
  min_idx = argmin(iterations_list[valid_indices])
  valid_omegas = omega_values[valid_indices]
  optimal_omega = valid_omegas[min_idx]
  min_iterations = minimum(iterations_list[valid_indices])

  println("Optimalna omega: ", optimal_omega)
  println("Minimalno število iteracij: ", min_iterations)
end

# Narišemo graf
omega_plot = plot(omega_values, iterations_list,
  xlabel="ω",
  ylabel="Število iteracij",
  title="Odvisnost hitrosti konvergence SOR metode od ω",
  marker=:circle,
  linewidth=2)

# Označimo optimalno točko
if any(valid_indices)
  scatter!(omega_plot, [optimal_omega], [min_iterations],
    color=:red, markersize=8, label="Optimalna ω")
end

savefig(omega_plot, "docs/src/assets/omega_convergence.png")
println("Graf shranjen v 'docs/src/assets/omega_convergence.png'")
