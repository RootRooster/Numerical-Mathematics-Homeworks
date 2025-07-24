using Homework1
using LinearAlgebra
using Test

@testset "Homework1.jl" begin

  @testset "Constructor Tests" begin
    A_dense = [
      10 0 -2;
      3 14 0;
      0 0 5
    ]
    sparse_A = RedkaMatrika(A_dense)

    @test length(sparse_A.V) == 3
    @test typeof(sparse_A) == RedkaMatrika{Int}
    @test sparse_A.V[1] == [10, -2]
    @test sparse_A.I[1] == [1, 3]

    @testset "Failure: Non-Square Matrix" begin
      A_nonsquare = [1 2 3; 4 5 6]
      @test_throws DimensionMismatch RedkaMatrika(A_nonsquare)
    end
  end

  @testset "Interface Tests" begin
    A_dense_float = [
      10.0 0.0 -2.0;
      3.0 14.0 0.0;
      0.0 0.0 5.0
    ]

    @testset "size" begin
      A_sparse = RedkaMatrika(A_dense_float)
      @test size(A_sparse) == (3, 3)
      # Edge case: empty matrix
      A_empty = RedkaMatrika(zeros(0, 0))
      @test size(A_empty) == (0, 0)
    end

    @testset "getindex (A[i, j])" begin
      A_sparse = RedkaMatrika(A_dense_float)
      # Reading an existing non-zero element
      @test A_sparse[1, 1] == 10.0
      @test A_sparse[3, 3] == 5.0

      # Reading a zero element (which is not explicitly stored)
      @test A_sparse[1, 2] == 0.0

      # Your existing bounds check test is already perfect
      @testset "Bounds Checking" begin
        @test_throws BoundsError A_sparse[4, 1]
        @test_throws BoundsError A_sparse[1, 4]
      end
    end

    @testset "setindex! (A[i, j] = value)" begin

      @testset "Update: Non-Zero to Non-Zero" begin
        A = RedkaMatrika(A_dense_float)
        A[1, 1] = 99.0
        @test A[1, 1] == 99.0
        @test A.V[1] == [99.0, -2.0] # Check internal state
      end

      @testset "Insertion: Zero to Non-Zero" begin
        A = RedkaMatrika(A_dense_float)
        A[1, 2] = 7.7 # Was zero, now non-zero
        @test A[1, 2] == 7.7
        @test A.V[1] == [10.0, 7.7, -2.0] # Check sorted insertion
        @test A.I[1] == [1, 2, 3]
      end

      @testset "Deletion: Non-Zero to Zero" begin
        A = RedkaMatrika(A_dense_float)
        A[1, 1] = 0.0
        @test A[1, 1] == 0.0
      end

      @testset "No-Op: Zero to Zero" begin
        A = RedkaMatrika(A_dense_float)
        V_before = deepcopy(A.V[1]) # Save state before operation

        A[1, 2] = 0.0 # This element was already zero

        @test A[1, 2] == 0.0
        @test A.V[1] == V_before # Ensure storage did not change
      end

      # Your bounds check is great, just giving it a clear spot
      @testset "Bounds Checking" begin
        A = RedkaMatrika(A_dense_float)
        @test_throws BoundsError (A[4, 1] = 100.0)
      end
    end

    @testset "Matrix-Vector Multiplication (*)" begin
      A_sparse = RedkaMatrika(A_dense_float)
      v = [1.0, 2.0, 3.0]

      # Test success by comparing to dense matrix result
      @test isapprox(A_sparse * v, A_dense_float * v)

      @testset "Failure: Dimension Mismatch" begin
        v_wrong_dim = [1.0, 2.0]
        @test_throws DimensionMismatch A_sparse * v_wrong_dim
      end
    end
  end

  @testset "SOR Method" begin
    # Setup for a known problem
    A_dense = [
      4.0 -1.0 -1.0 0.0;
      -1.0 4.0 0.0 -1.0;
      -1.0 0.0 4.0 -1.0;
      0.0 -1.0 -1.0 4.0
    ]
    b = [2.0, 2.0, 2.0, 2.0]
    x0 = zeros(4)
    # The exact solution to this system is a vector of ones.
    x_exact = [1.0, 1.0, 1.0, 1.0]
    omega = 1.2

    A_sparse = RedkaMatrika(A_dense)

    @testset "Successful Convergence" begin
      x, it = sor(A_sparse, b, x0, omega, tol=1e-12)

      @test it > 0
      @test it < 100
      @test isapprox(x, x_exact, atol=1e-10)
    end

    @testset "Input Validation and Errors" begin
      b_wrong = [1.0, 2.0]
      x0_wrong = [1.0, 2.0, 3.0]
      @test_throws DimensionMismatch sor(A_sparse, b_wrong, x0, omega)
      @test_throws DimensionMismatch sor(A_sparse, b, x0_wrong, omega)

      A_zero_diag_dense = [4.0 -1.0; -1.0 0.0] # Zero on diagonal
      A_zero_diag = RedkaMatrika(A_zero_diag_dense)
      b2 = [1.0, 1.0]
      x02 = [0.0, 0.0]
      @test_throws ErrorException sor(A_zero_diag, b2, x02, omega)
    end

    @testset "Parameter Effects (tol)" begin
      # A loose tolerance should result in fewer iterations.
      _, it_loose = sor(A_sparse, b, x0, omega, tol=1e-2)

      # A tight tolerance should result in more iterations.
      _, it_tight = sor(A_sparse, b, x0, omega, tol=1e-10)

      @test it_loose > 0
      @test it_tight > it_loose
    end

    @testset "Convergence Failure" begin
      # This matrix is not diagonally dominant and will likely not converge.
      A_non_conv_dense = [1.0 2.0; 3.0 1.0]
      A_non_conv = RedkaMatrika(A_non_conv_dense)
      b2 = [1.0, 1.0]
      x02 = [0.0, 0.0]

      @test_warn "SOR method did not converge" begin
        # We set a very low max_iter to force the failure.
        x, it = sor(A_non_conv, b2, x02, 1.1, max_iter=5)

        @test it == 5
        x_true_non_conv = [1 / 5, 2 / 5]
        @test !isapprox(x, x_true_non_conv)
      end
    end
  end

end
