module Homework1
using LinearAlgebra

export RedkaMatrika, sor

"""
    RedkaMatrika{T<:Number}

A custom sparse matrix data type.
"""
struct RedkaMatrika{T<:Number}
  V::Vector{Vector{T}}   # Stores the non-zero values
  I::Vector{Vector{Int}} # Stores the column indices
end

"""
    RedkaMatrika(A::Matrix{T}) where {T<:Number}

Constructor for `RedkaMatrika` from a dense matrix `A`.
"""
function RedkaMatrika(A::Matrix{T}) where {T<:Number}
  rows, cols = size(A)
  if rows != cols
    throw(DimensionMismatch("Input matrix must be square."))
  end
  n = rows

  V_vecs = [Vector{T}() for _ in 1:n]
  I_vecs = [Vector{Int}() for _ in 1:n]

  for i in 1:n
    for j in 1:n
      if A[i, j] != zero(T)
        push!(V_vecs[i], A[i, j])
        push!(I_vecs[i], j)
      end
    end
  end

  return RedkaMatrika(V_vecs, I_vecs)
end


import Base: getindex, setindex!, size, *

"""
    size(A::RedkaMatrika) -> Tuple{Int, Int}

Returns the dimensions of the sparse matrix as a tuple (n, n).
"""
size(A::RedkaMatrika) = (length(A.V), length(A.V))

"""
    getindex(A::RedkaMatrika{T}, i::Int, j::Int) where T

Returns the value of the element at position (i, j). Enables `A[i, j]` syntax.
"""
function getindex(A::RedkaMatrika{T}, i::Int, j::Int) where {T<:Number}
  @boundscheck checkbounds(1:length(A.V), i)
  @boundscheck checkbounds(1:length(A.V), j)

  # Find the idx value 
  row_indices = A.I[i]
  idx = searchsortedfirst(row_indices, j)

  if idx <= length(row_indices) && row_indices[idx] == j
    # Element is a stored non-zero value
    return A.V[i][idx]
  else
    # Element is not stored, so it's zero
    return zero(T)
  end
end

"""
    setindex!(A::RedkaMatrika{T}, value::T, i::Int, j::Int) where T

Sets the value of the element at position (i, j). Enables `A[i, j] = value` syntax.
"""
function setindex!(A::RedkaMatrika{T}, value::T, i::Int, j::Int) where {T<:Number}
  @boundscheck checkbounds(1:length(A.V), i)
  @boundscheck checkbounds(1:length(A.V), j)

  # Find the column index `j` in the i-th row
  idx_in_row = findfirst(isequal(j), A.I[i])

  if idx_in_row === nothing
    # Element (i, j) is not currently stored (it's zero).
    if !iszero(value)
      # Insert the new non zero element at the sorted position.
      insert_pos = searchsortedfirst(A.I[i], j)
      insert!(A.I[i], insert_pos, j)
      insert!(A.V[i], insert_pos, value)
    end
    # If the value is zero, do nothing.
  else
    # Element (i, j) already exists in the sparse structure.
    if iszero(value)
      # The value is being set to zero, so remove it from storage.
      deleteat!(A.I[i], idx_in_row)
      deleteat!(A.V[i], idx_in_row)
    else
      # Update the existing non-zero value.
      A.V[i][idx_in_row] = value
    end
  end
  return A
end


"""
    *(A::RedkaMatrika{T}, v::Vector{T}) where T -> Vector{T}

Computes the product of a sparse matrix `A` and a vector `v`.
"""
function *(A::RedkaMatrika{T}, v::Vector{T}) where {T<:Number}
  n = length(A.V)
  if n != length(v)
    throw(DimensionMismatch("Matrix dimensions $(size(A)) do not match vector length $(length(v))"))
  end

  # Initialize the result vector with zeros
  result = zeros(T, n)

  # Iterate only over the non-zero elements
  for i in 1:n # Iterate over the matrix rows
    # Sum for the i-th row of the result
    row_sum = zero(T)
    for k in 1:length(A.I[i])
      col_idx = A.I[i][k]
      value = A.V[i][k]
      row_sum += value * v[col_idx]
    end
    result[i] = row_sum
  end

  return result
end

"""
    sor(A, b, x0, omega; tol=1e-10, max_iter=1000) -> (Vector, Int)

Solves the linear system Ax = b using the Successive Over-Relaxation (SOR) iterative method.

# Arguments
- `A::RedkaMatrika`: The square, sparse matrix of the system.
- `b::Vector`: The right-hand side vector.
- `x0::Vector`: The initial guess for the solution.
- `omega::Real`: The relaxation parameter (typically 0 < omega < 2).
- `tol::Real=1e-10`: The tolerance for stopping; iteration stops when the infinity norm of the residual is less than `tol`.
- `max_iter::Int=1000`: The maximum number of iterations.

# Returns
- `Vector`: The computed approximate solution `x`.
- `Int`: The number of iterations performed.
"""
function sor(A::RedkaMatrika{T}, b::Vector{T}, x0::Vector{T}, omega::Real; tol::Real=1e-10, max_iter::Int=1000) where {T<:Number}
  n = length(A.V)
  if n != length(b) || n != length(x0)
    throw(DimensionMismatch("Dimenzije matrike A ter vektorjev b in x0 se ne ujemajo."))
  end

  # Set the initial x
  x = copy(x0)

  # Main SOR iteration loop
  for it in 1:max_iter

    # Update each component x[i] of the solution vector
    for i in 1:n
      sigma = zero(T)
      diag_val = zero(T)
      found_diag = false

      # Iterate over non-zero elements of row 'i'
      for idx in 1:length(A.I[i])
        j = A.I[i][idx]
        val = A.V[i][idx]

        # Separate the digonal values from the off-diagonal sum (sigma)
        if i == j
          diag_val = val
          found_diag = true
        else
          sigma += val * x[j]
        end
      end

      # SOR requres non-zero diagonal elements for the division
      if !found_diag || iszero(diag_val)
        error("Matrix has missing or zero diagonal element in row $i. SOR method is NOT applicable.")
      end

      # Store old x[i] value 
      x_i_old = x[i]

      # Apply core SOR iteration formula
      term = (b[i] - sigma) / diag_val
      x[i] = (1 - omega) * x_i_old + omega * term
    end

    # After full swap (all x[i] are updated), check the stopping condition
    residual = A * x - b
    if norm(residual, Inf) < tol
      return x, it # Success: converged within tolerance
    end
  end
  @warn("SOR method did not converge within $max_iter iterations")
  return x, max_iter
end

end # End of module Homework1
