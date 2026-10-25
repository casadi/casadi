# native_inputs.jl -- Julia numbers, arrays and ranges are accepted wherever a DM
# is, and Opti values come back as native Julia values.

@testset "native inputs: Function calls" begin
    x = sym(MX, "x", 2); p = sym(MX, "p")
    F = C.Function("F", [x, p], [x * p])
    @test Matrix(F([1, 2], 3)) == reshape([3.0, 6.0], 2, 1)
    @test Matrix(F(1:2, 3)) == reshape([3.0, 6.0], 2, 1)
    @test Matrix(F(range(1, 2; length = 2), [3])) == reshape([3.0, 6.0], 2, 1)
    @test F(x, [3]) isa MX
    @test Matrix(F.map(2)([1 3; 2 4], [1 2])) == [1.0 6.0; 2.0 8.0]
    E = C.Function("E", [sym(MX, "e", 0, 0)], [MX(1.0)])
    @test Float64(E([])) == 1.0
end

@testset "native inputs: conversions" begin
    @test Matrix(DM(range(-1, 1; length = 3))) == reshape([-1.0, 0.0, 1.0], 3, 1)
    @test Matrix(DM([1 2; 3 4])) == [1.0 2.0; 3.0 4.0]
    @test size(DM([1.0, 2.0]')) == (1, 2)
    @test SX(1:3) isa SX
    @test MX(0:2) isa MX
    # any Real converts like a Float64
    @test Float64(DM(pi)) == Float64(pi)
    @test Float64(DM(1.5f0)) == 1.5
    @test Float64(DM(1 // 4)) == 0.25
    @test Float64(DM(big(2.5))) == 2.5
    x = sym(SX, "x")
    @test Float64(C.Function("f", [x], [2 * x])(pi)) == 2 * Float64(pi)
end

@testset "native inputs: free functions keep the casadi type" begin
    x = sym(MX, "x", 2); s = sym(SX, "s", 2)
    @test C.fmax(x, 0) isa MX
    @test C.fmax(s, 0) isa SX
    @test Matrix(C.fmax([1, -2], 0)) == reshape([1.0, 0.0], 2, 1)
    @test C.fmax(1, 2) == 2.0
    @test C.if_else(x[1] >= 0, 1 / 0.9, 0.8) isa MX
    @test C.logsumexp(x, 1) isa MX
    @test C.mtimes([1.0 2.0; 3.0 4.0], s) isa SX
end

@testset "native inputs: arithmetic and comparisons" begin
    x = sym(MX, "x", 2); s = sym(SX, "s", 2)
    @test (x - [1.0, 2.0]) isa MX
    @test ([1, 2] + x) isa MX
    @test ([1.0 2.0; 3.0 4.0] * s) isa SX
    @test ([1.0 2.0; 3.0 4.0] \ s) isa SX
    @test Matrix(DM([1, 2]) + [1, 1]) == reshape([2.0, 3.0], 2, 1)
    @test (x == [1, 2]) isa MX
    @test (s >= zeros(2)) isa SX
    @test dot([1, 2], x) isa MX
end

@testset "native inputs: concatenation, casadi operand anywhere" begin
    v = sym(MX, "v"); s = sym(SX, "s"); x = sym(MX, "x", 2)
    @test size([1 2 3 4 5 6 7 8 9 10 11 v]) == (1, 12)
    @test size([[1.0 2.0; 3.0 4.0] x; x' 0]) == (3, 3)
    w = [sym(MX, "w$i") for i in 1:300]
    @test size(vcat(w...)) == (300, 1)
    @test size(hcat(w...)) == (1, 300)
    @test size(reduce(vcat, w)) == (300, 1)
    @test size([w...;]) == (300, 1)
    @test [Any[1]; 2] isa Vector{Any}       # no casadi operand: still Base
    @test [0; 0; v] isa MX
    @test size([0; 0; v]) == (3, 1)
    @test size([0 0 v]) == (1, 3)
    @test size([1 v; 0 0]) == (2, 2)
    @test [1; 2; 3; 4; 5; 6; 7; 8; s] isa SX
    @test size([v; [1, 2]]) == (3, 1)
    @test [0; 0; DM(1)] isa DM
    @test [1; 2] isa Vector{Int}            # all-Julia cats stay Base
    @test size([DM([1.0, 2.0]) [3.0, 4.0]]) == (2, 2)
    @test size(hcat(DM[1, 2]...)) == (1, 2)  # a Julia vector of casadi values
    @test [v; [1, 2]] isa MX
end

HAS_IPOPT && @testset "native inputs: Opti" begin
    opti = Opti()
    y = variable(opti, 3); p = parameter(opti)
    opti.set_value(p, 2)
    opti.set_initial(y, range(-1, 1; length = 3))
    opti.subject_to(y[1:2] == [1, 2])
    minimize(opti, C.sumsqr(y - p))
    solver(opti, "ipopt", Dict("print_time" => false, "ipopt" => Dict("print_level" => 0)))
    sol = solve(opti)
    @test sol.value(y) isa Vector{Float64}
    @test sol.value(y) ≈ [1.0, 2.0, 2.0] atol = 1e-6
    @test sol.value(y[3]) isa Float64
    @test sol.value(y') isa Vector{Float64}
    @test sol.value([y y]) isa Matrix{Float64}
    @test value(sol, y) ≈ sol.value(y)
    @test opti.debug.value(y) ≈ sol.value(y)
    # Opti members Python exposes as properties are values here too, not bound methods
    @test opti.f isa MX && opti.x isa MX && opti.g isa MX && opti.lam_g isa MX
    @test opti.nx == 3 && opti.ng == 2
    @test sol.opti isa C.Opti
    @test sol.value(C.gradient(opti.f + opti.lam_g' * opti.g, opti.x)) ≈ zeros(3) atol = 1e-6
    @test opti.to_function("G", [p], [y]) isa C.Function    # real methods stay methods
end

@testset "tuple dims, untyped empty vectors, typemax(Int)" begin
    @test size(C.dense(Sparsity, (2, 3))) == (2, 3)
    @test size(C.nan(DM, (2, 1))) == (2, 1)
    @test Matrix(C.repmat(DM(7), (1, 3))) == [7 7 7]
    @test C.nnz(C.rowcol(Sparsity, [], [1], 3, 4)) == 0
    @test size(C.triplet(DM, [], [], [], 2, 2)) == (2, 2)
    @test C.uni_coloring(C.diag(Sparsity, 3), Sparsity(), typemax(Int)) isa Sparsity
end
