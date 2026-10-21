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
