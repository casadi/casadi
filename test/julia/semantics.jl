# semantics.jl -- Julia's equality, comparison, zero/one and indexing conventions on casadi values.

@testset "equality and comparison" begin
    x = sym(SX, "x"); v = sym(SX, "v", 3); M = sym(MX, "M", 2)
    D = DM([1 2 3; 4 5 6])
    # SX/MX: == and != build constraints; isequal/hash are structural identity
    @test (x != v[1]) isa SX
    @test isequal(x, x) && !isequal(x, sym(SX, "x"))
    @test isequal(M, M) && hash(M) == hash(M)
    @test Dict(x => 1, v => 2)[v] == 2
    @test length(unique([x, x, v[1]])) == 2
    @test iszero(SX(0)) && !iszero(x) && iszero(MX(0))
    # DM compares elementwise into a DM, like SX/MX; isequal/hash/iszero are Bool
    @test Matrix(D == [1 2 3; 4 5 6]) == ones(2, 3)
    @test (D <= 3) isa DM && Matrix(D <= 3) == [1.0 1.0 1.0; 0.0 0.0 0.0]
    @test Matrix(2 < D) == [0.0 0.0 1.0; 1.0 1.0 1.0]
    @test Float64(DM(1) != DM(2)) == 1.0 && Float64(DM(3) >= DM(3)) == 1.0
    @test isequal(DM(2), DM(2)) && !isequal(D, D') && length(Set([DM(1), DM(1), DM(2)])) == 2
    @test iszero(DM(0)) && !iszero(D)
    @test Matrix(D .> 2) == Matrix(D > 2)
    @test Matrix(D .== 2) == [0.0 1.0 0.0; 0.0 0.0 0.0]
    @test (x .< v) isa SX && size(x .< v) == (3, 1)
end

@testset "zero, one, Vector, ifelse" begin
    A = sym(SX, "A", 2, 2)
    @test size(zero(A)) == (2, 2) && iszero(zero(A))
    @test C.is_equal(one(A), C.eye(SX, 2))
    @test Matrix(one(DM([1 2; 3 4]))) == [1.0 0.0; 0.0 1.0]
    @test_throws DimensionMismatch one(DM([1 2 3]))
    v = sym(MX, "v", 3)
    @test Vector(v) isa Vector{MX} && length(Vector(v')) == 3
    @test_throws DimensionMismatch Vector(sym(SX, "B", 2, 2))
    x = sym(SX, "x")
    f = C.Function("f", [x], [ifelse(x > 0, x, -x)])
    @test Float64(f(-3)) == 3.0
    @test Matrix(ifelse(DM([1, 0]), 5, 7)) == reshape([5.0, 7.0], 2, 1)
end

@testset "indexing with integer vectors and Bool masks" begin
    D = DM([1 2 3; 4 5 6])
    @test Matrix(D[[2, 1]]) == reshape([4.0, 1.0], 2, 1)
    @test Matrix(D[[true, false, true, false, false, true]]) == reshape([1.0, 2.0, 6.0], 3, 1)
    @test Matrix(D[Bool[1 0 1; 0 1 0]]) == reshape([1.0, 5.0, 3.0], 3, 1)
    @test Matrix(D[1, [3, 1]]) == [3.0 1.0]
    @test Matrix(D[[2, 1], 2:3]) == [5.0 6.0; 2.0 3.0]
    @test Matrix(D[:, [true, false, true]]) == [1.0 3.0; 4.0 6.0]
    @test_throws BoundsError D[[7]]
    @test_throws BoundsError D[[0]]
    @test_throws BoundsError D[[true, false]]
    v = sym(SX, "v", 4)
    @test C.is_equal(v[[4, 1]], [v[4]; v[1]])
    E = DM([1 2 3; 4 5 6]); E[[1, 6]] = [10, 60]
    @test Matrix(E) == [10.0 2.0 3.0; 4.0 5.0 60.0]
    E = DM([1 2 3; 4 5 6]); E[1, [3, 1]] = [30 10]     # mixed Slice/vector set
    @test Matrix(E) == [10.0 2.0 30.0; 4.0 5.0 6.0]
    E = DM([1 2 3; 4 5 6]); E[E .> 4] = 0                # a DM from .> is a mask
    @test Matrix(E) == [1.0 2.0 3.0; 4.0 0.0 0.0]
    @test Matrix(D[D .>= 5]) == reshape([5.0, 6.0], 2, 1)
    @test Matrix(D[D >= 5]) == reshape([5.0, 6.0], 2, 1)
end

@testset "no silent departures from Julia's array semantics" begin
    A = DM([1.0 2; 3 4]); B = DM([2.0 0; 1 3])
    a = [1.0 2; 3 4]; b = [2.0 0; 1 3]
    S = sym(SX, "S", 2, 2); v = sym(MX, "v", 2); x = sym(SX, "x")
    # / and ^ are the matrix operations; the dotted forms are elementwise
    @test Matrix(A / B) ≈ a / b && Matrix(A ./ B) ≈ a ./ b
    @test Matrix(A^2) ≈ a^2 && Matrix(A .^ 2) ≈ a .^ 2 && Matrix(A^-1) ≈ inv(a)
    @test C.nnz(C.eye(DM, 3)^3) == 3                                 # mpower keeps the sparsity
    @test C.is_equal(C.simplify(S^2 - S * S), SX.zeros(2, 2))
    @test_throws ArgumentError A^0.5                                 # Julia's would be a matrix root
    @test_throws ArgumentError 2^A                                   # ... and a matrix exponential
    @test_throws ArgumentError S^0.5
    @test_throws DimensionMismatch v^2
    @test_throws ArgumentError 2 / A
    @test Float64(DM(6) / 2) == 3.0 && (x / 2) isa SX && (x^2) isa SX && Float64(DM(2)^0.5) ≈ sqrt(2)
    # matrix functions: casadi has none, so a non-scalar throws instead of going elementwise
    for M in (A, S, v)
        @test_throws ArgumentError sin(M)
        @test_throws ArgumentError exp(M)
    end
    @test sin(x) isa SX && size(sin.(S)) == (2, 2) && size(exp.(v)) == (2, 1)
    @test Matrix(exp.(A)) ≈ exp.(a)
    @test Matrix(abs(DM([-1, 2]))) == reshape([1.0, 2.0], 2, 1)      # no array meaning: elementwise
    # max/min of arrays are Julia's lexicographic pick: scalars only, elementwise with dots
    @test_throws ArgumentError max(A, B)
    @test_throws ArgumentError max(DM([1.0, 2]), DM([3.0, 0]))
    @test_throws ArgumentError max(S, 0)
    @test Matrix(max.(A, B)) == max.(a, b) && Float64(max(DM(1), 2)) == 2.0
    # an arbitrary f. is applied per element
    f(z) = z^2 + 1
    @test Matrix(f.(A)) ≈ f.(a)
    @test C.is_equal(C.simplify(f.(S) - (S .^ 2 .+ 1)), SX.zeros(2, 2))
    # shapes follow casadi: A[i, :] is a row, a linear range keeps a vector's orientation
    r = DM([1.0 2 3])
    @test size(A[1, :]) == (1, 2) && size(r[1:2]) == (1, 2) && size(r[:]) == (1, 3)
    @test size(r[[1, 3]]) == (1, 2) && size(A[2:3]) == (2, 1)
    R = DM([1.0 2 3]); R[1:2] = [10, 20]
    @test Matrix(R) == [10.0 20.0 3.0]
    E = DM([1.0 2; 3 4]); E[2, :] = [7, 8]
    @test Matrix(E) == [1.0 2.0; 7.0 8.0]
end
