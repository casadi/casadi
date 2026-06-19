# multiplication.jl -- capability port of test/javascript/multiplication.js +
# test/python/multiplication.py subset.  mtimes (matrix product) vs elementwise
# times, matrix*vector, transpose, kron, dot, cross, outer (rank-1),
# associativity, identity, sparse*dense, mixed DM/SX/MX promotion in products,
# and A \ b consistency with inv(A)*b.  All assertions are analytic.

# ---- mtimes vs elementwise on DM ----
@testset "mtimes vs elementwise times" begin
    A = DM([1.0 2.0 3.0; 4.0 5.0 6.0])     # 2x3
    B = DM([7.0 8.0; 9.0 10.0; 11.0 12.0]) # 3x2
    C2 = C.mtimes(A, B)                      # 2x2
    @test (C.size1(C2), C.size2(C2)) == (2, 2)
    @test Matrix{Float64}(C2) == [58.0 64.0; 139.0 154.0]
    # C.times is the explicit elementwise (Hadamard) product.
    P = C.times(DM([1.0 2.0; 3.0 4.0]), DM([5.0 6.0; 7.0 8.0]))
    @test Matrix{Float64}(P) == [5.0 12.0; 21.0 32.0]
    @test Matrix{Float64}(C.times(DM([2.0 3.0]), DM([4.0 5.0]))) == [8.0 15.0]
    # broadcasting `.*` is elementwise (Hadamard), distinct from `*`/mtimes.
    @test Matrix{Float64}(DM([1.0 2.0; 3.0 4.0]) .* DM([5.0 6.0; 7.0 8.0])) == [5.0 12.0; 21.0 32.0]
    @test Matrix{Float64}(DM([1.0 2.0; 3.0 4.0]) * DM([5.0 6.0; 7.0 8.0])) == [19.0 22.0; 43.0 50.0]
    # mtimes != times: for 2x2 the shapes match but values differ.
    M22 = DM([1.0 2.0; 3.0 4.0]); N22 = DM([5.0 6.0; 7.0 8.0])
    @test Matrix{Float64}(C.mtimes(M22, N22)) == [19.0 22.0; 43.0 50.0]
    @test Matrix{Float64}(C.times(M22, N22)) == [5.0 12.0; 21.0 32.0]
end

# ---- matrix * vector ----
@testset "matrix * vector" begin
    A = DM([1.0 2.0; 3.0 4.0]); v = DM([5.0, 6.0])
    @test Vector{Float64}(C.mtimes(A, v)) == [17.0, 39.0]   # [1*5+2*6, 3*5+4*6]
    # row-vector * matrix.
    r = DM([1.0 2.0])
    @test Matrix{Float64}(C.mtimes(r, A)) == [7.0 10.0]     # [1+6, 2+8]
end

# ---- transpose interplay: (A B)' = B' A' ----
@testset "transpose of a product" begin
    A = DM([1.0 2.0 3.0; 4.0 5.0 6.0]); B = DM([7.0 8.0; 9.0 10.0; 11.0 12.0])
    lhs = C.transpose(C.mtimes(A, B))
    rhs = C.mtimes(C.transpose(B), C.transpose(A))
    @test Matrix{Float64}(lhs) == Matrix{Float64}(rhs)
    @test (C.size1(lhs), C.size2(lhs)) == (2, 2)
end

# ---- kron (Kronecker product) ----
@testset "kron" begin
    K = C.kron(DM([1.0 0.0; 0.0 1.0]), DM([1.0 2.0; 3.0 4.0]))
    @test (C.size1(K), C.size2(K)) == (4, 4)
    @test Matrix{Float64}(K) == [1.0 2.0 0.0 0.0; 3.0 4.0 0.0 0.0;
                                 0.0 0.0 1.0 2.0; 0.0 0.0 3.0 4.0]
    # kron(a,b) shape = (rows_a*rows_b, cols_a*cols_b).
    K2 = C.kron(DM([1.0 2.0]), DM([3.0, 4.0]))   # (1x2) kron (2x1) -> 2x2
    @test (C.size1(K2), C.size2(K2)) == (2, 2)
    @test Matrix{Float64}(K2) == [3.0 6.0; 4.0 8.0]
end

# ---- dot / cross / outer (rank-1) ----
@testset "dot / cross / outer" begin
    a = DM([1.0, 2.0, 3.0]); b = DM([4.0, 5.0, 6.0])
    @test Float64(C.dot(a, b)) == 32.0
    # cross product of unit axes: e1 x e2 = e3.
    @test Vector{Float64}(C.cross(DM([1.0, 0.0, 0.0]), DM([0.0, 1.0, 0.0]))) == [0.0, 0.0, 1.0]
    @test Vector{Float64}(C.cross(a, b)) == [-3.0, 6.0, -3.0]  # [2*6-3*5, 3*4-1*6, 1*5-2*4]
    # outer product a b' is rank-1: column*row via mtimes.
    O = C.mtimes(a, C.transpose(b))
    @test (C.size1(O), C.size2(O)) == (3, 3)
    @test Matrix{Float64}(O) == [4.0 5.0 6.0; 8.0 10.0 12.0; 12.0 15.0 18.0]
end

# ---- associativity and identity ----
@testset "associativity & identity" begin
    A = DM([1.0 2.0; 0.0 1.0]); B = DM([1.0 0.0; 3.0 1.0]); D = DM([2.0 1.0; 1.0 2.0])
    left = C.mtimes(C.mtimes(A, B), D)
    right = C.mtimes(A, C.mtimes(B, D))
    @test Matrix{Float64}(left) == Matrix{Float64}(right)     # (AB)D == A(BD)
    I2 = C.eye(C.DM, 2)
    @test Matrix{Float64}(C.mtimes(A, I2)) == Matrix{Float64}(A)
    @test Matrix{Float64}(C.mtimes(I2, A)) == Matrix{Float64}(A)
end

# ---- sparse * dense yields the correct dense result ----
@testset "sparse * dense" begin
    A = C.sparsify(DM([1.0 0.0 0.0; 1.0 1.0 0.0; 1.0 1.0 1.0]))  # lower-triangular
    @test C.nnz(C.sparsity(A)) == 6                               # 6 stored nonzeros
    B = DM([1.0 2.0; 3.0 4.0; 5.0 6.0])
    Cd = C.mtimes(A, B)
    @test (C.size1(Cd), C.size2(Cd)) == (3, 2)
    @test Matrix{Float64}(Cd) == [1.0 2.0; 4.0 6.0; 9.0 12.0]
end

# ---- mixed DM/SX/MX promotion in products ----
@testset "mixed-type promotion in mtimes" begin
    # SX symbolic mtimes evaluated numerically.
    As = sym(SX, "A", 2, 2); xs = sym(SX, "x", 2)
    fsx = C.Function("fsx", [As, xs], [C.mtimes(As, xs)])
    rsx = fsx(DM([1.0 2.0; 3.0 4.0]), DM([5.0, 6.0]))
    @test Vector{Float64}(rsx) == [17.0, 39.0]
    # MX symbolic mtimes evaluated numerically (mirrors multiplication.js test_mtimes_mx_eval).
    Am = sym(MX, "A", 2, 3); Bm = sym(MX, "B", 3, 2)
    fmx = C.Function("fmx", [Am, Bm], [C.mtimes(Am, Bm)])
    rmx = fmx(DM([1.0 2.0 3.0; 4.0 5.0 6.0]), DM([7.0 8.0; 9.0 10.0; 11.0 12.0]))
    @test Matrix{Float64}(rmx) == [58.0 64.0; 139.0 154.0]
    # operator-* promotes DM into the SX product.
    y = sym(SX, "y")
    prod = DM([2.0]) * y                                          # DM * SX -> SX
    @test prod isa SX
    @test abs(Float64((C.Function("g", [y], [prod]))(DM(3.0))) - 6.0) < 1e-12
end

# ---- A \ b consistency: solve(A,b) == inv(A) * b ----
@testset "solve vs inv*b consistency" begin
    A = DM([3.0 1.0; 1.0 2.0]); b = DM([9.0, 8.0])
    x_solve = C.solve(A, b, "qr", Dict{String,Any}())
    x_inv = C.mtimes(C.inv(A), b)
    @test maximum(abs.(Vector{Float64}(x_solve) .- Vector{Float64}(x_inv))) < 1e-10
    # verify the residual A x - b is ~0.
    resid = Vector{Float64}(C.minus(C.mtimes(A, x_solve), b))
    @test maximum(abs.(resid)) < 1e-10
    # analytic: A x = b -> x = [2, 3] (3*2+1*3=9, 1*2+2*3=8).
    @test maximum(abs.(Vector{Float64}(x_solve) .- [2.0, 3.0])) < 1e-10
end

# ---- is_compactible (mirror multiplication.js test_is_compactible_dense) ----
@testset "dense sparsity compactness" begin
    sp = C.dense(C.Sparsity, 3, 4)
    @test C.is_dense(sp)
    @test C.nnz(sp) == 12
end
