# base_math.jl -- Julia's own reductions and math functions work on casadi values,
# with Julia's semantics, and keep the casadi type.
using LinearAlgebra: norm, det, tr

@testset "Base/LinearAlgebra math on casadi values" begin
    x = sym(MX, "x", 3); s = sym(SX, "s", 2); A = DM([1 2 3; 4 5 6])
    @test sum(x) isa MX && sum(s) isa SX
    @test Float64(sum(A)) == 21.0                          # all elements, as for arrays
    @test Matrix(sum(A; dims = 1)) == [5.0 7.0 9.0]
    @test Matrix(sum(A; dims = 2)) == reshape([6.0, 15.0], 2, 1)
    @test diff(x) isa MX && size(diff(x)) == (2, 1)
    @test Matrix(diff(A; dims = 2)) == [1.0 1.0; 1.0 1.0]
    @test Matrix(diff(DM([1 4 9]))) == [3.0 5.0]           # a row diffs along the row
    @test_throws ArgumentError diff(A)                     # a matrix needs dims, as in Julia
    @test Matrix(cumsum(A; dims = 2)) == [1.0 3.0 6.0; 4.0 9.0 15.0]
    @test Float64(maximum(A)) == 6.0 && Float64(minimum(A)) == 1.0
    @test max.(x, 0) isa MX && max.(0, x) isa MX && min.(s, [1, 2]) isa SX
    @test max(x[1], 0) isa MX && min(0, s[2]) isa SX
    @test max(x[1], x[2], 0) isa MX
    @test max(1, 2) == 2                                    # plain numbers stay Base
    @test norm(x) isa MX && norm(x, 1) isa MX && norm(x, Inf) isa MX
    @test Float64(norm(DM([3, 4]))) ≈ 5.0
    @test Float64(norm(A)) ≈ sqrt(91)                       # Frobenius for a matrix
    @test norm([3, 4]) == 5.0
    M = DM([1 2; 3 4])
    @test Matrix(inv(M)) ≈ [-2.0 1.0; 1.5 -0.5]
    @test Float64(det(M)) ≈ -2.0 && Float64(tr(M)) == 5.0
    @test atan(x[1], 1) isa MX && hypot(x[1], x[2]) isa MX
    @test Float64(log10(DM(100))) ≈ 2.0
end

@testset "logic: & | ! on scalars, .& .| .! elementwise" begin
    x = sym(SX, "x"); y = sym(SX, "y")
    @test C.str(x & y) == "(x&&y)" && C.str(x | y) == "(x||y)" && C.str(!x) == "(!x)"
    @test Bool((DM(2) > 1) & (DM(1) > 3)) === false
    @test Matrix(DM([1, 0]) .& DM([1, 1])) == reshape([1.0, 0.0], 2, 1)
    @test Matrix(DM([1, 0]) .| [0, 0]) == reshape([1.0, 0.0], 2, 1)
    @test Matrix(.!DM([1, 0])) == reshape([0.0, 1.0], 2, 1)
    @test_throws ArgumentError DM([1, 0]) & DM([1, 1])
    @test_throws ArgumentError !sym(MX, "z", 2)
    @test Matrix(C.logic_and(DM([1, 0]), DM([1, 1]))) == reshape([1.0, 0.0], 2, 1)
    @test Float64(C.logic_all(DM([1, 1]))) == 1
end

@testset "rem / mod / fld / cld / div / copysign with Julia's rounding" begin
    x = sym(SX, "x")
    f = C.Function("f", [x], [mod(x, 3), x % 3, rem(x, -3), fld(x, 3), cld(x, 3), x ÷ 3, copysign(x, -1), mod(x, -3)])
    for v in (-7.5, 7.5, -3.0, 4.0)
        @test [Float64(o) for o in f(v)] == [mod(v, 3), v % 3, rem(v, -3), fld(v, 3), cld(v, 3), v ÷ 3, copysign(v, -1), mod(v, -3)]
    end
    @test Matrix(mod.(DM([-1, 5]), 3)) == reshape([2.0, 2.0], 2, 1)
    @test_throws ArgumentError mod(DM([1, 2]), 3)
end

@testset "reshape with : and with a Sparsity" begin
    A = DM([1 2 3; 4 5 6])
    @test Matrix(reshape(A, :, 1)) == reshape([1 2 3; 4 5 6], :, 1)
    @test Matrix(reshape(A, 3, :)) == reshape([1 2 3; 4 5 6], 3, :)
    @test Matrix(reshape(A, (:, 2))) == reshape([1 2 3; 4 5 6], (:, 2))
    @test_throws DimensionMismatch reshape(A, :, 4)
    x = sym(SX, "x", 2, 2)
    @test size(reshape(x, C.dense(Sparsity, 4, 1))) == (4, 1)
end
