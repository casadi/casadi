# misc.jl -- capability port of test/javascript/misc.js + test/python/misc.py
# subset.  Meta/version strings, global-options round-trip, builder helpers
# (eye/diag/repmat/blockcat/veccat/linspace), type predicates (is_equal/is_op),
# string/printing, and linear-algebra value functions (dot/norms/trace/det/inv)
# with analytic checks.  The JS twin has no meta/global-options/linalg coverage.

# ---- meta strings: version / git_revision / compiler_id ----
@testset "CasadiMeta version strings" begin
    v = C.version(C.CasadiMeta)
    @test v isa AbstractString && occursin(r"^\d+\.\d+", v)   # "3.7.2+": major.minor
    rev = C.git_revision(C.CasadiMeta)
    @test rev isa AbstractString
    # a git sha: non-empty and purely hex (not just "length >= 7", which a
    # non-hex 7-char string would also pass); a full sha is 40 hex chars.
    @test !isempty(rev) && occursin(r"^[0-9a-fA-F]+$", rev) && length(rev) >= 7
    @test C.compiler_id(C.CasadiMeta) isa AbstractString
    @test !isempty(C.compiler_id(C.CasadiMeta))               # e.g. "GNU"
end

# ---- global options: getMaxNumDir / setMaxNumDir round-trip ----
@testset "GlobalOptions max-num-dir round-trip" begin
    orig = C.getMaxNumDir(C.GlobalOptions)
    @test orig isa Integer
    C.setMaxNumDir(C.GlobalOptions, 99)
    @test C.getMaxNumDir(C.GlobalOptions) == 99
    C.setMaxNumDir(C.GlobalOptions, orig)                     # restore
    @test C.getMaxNumDir(C.GlobalOptions) == orig
end

# ---- builders: eye / diag / repmat / blockcat / veccat / linspace ----
@testset "eye / ones / zeros (type-dispatched)" begin
    @test Matrix{Float64}(C.eye(C.DM, 3)) == [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test Matrix{Float64}(C.ones(C.DM, 2, 3)) == ones(2, 3)
    @test Matrix{Float64}(C.zeros(C.DM, 2, 3)) == zeros(2, 3)
    @test (C.size1(C.DM(C.eye(C.DM, 4))), C.size2(C.DM(C.eye(C.DM, 4)))) == (4, 4)
end

@testset "diag: vector<->matrix" begin
    @test Matrix{Float64}(C.diag(DM([1.0, 2.0, 3.0]))) ==
          [1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0]      # vector -> diagonal matrix
    @test Vector{Float64}(C.diag(DM([2.0 1.0; 1.0 3.0]))) == [2.0, 3.0]  # matrix -> diagonal
end

@testset "repmat / linspace" begin
    @test Matrix{Float64}(C.repmat(DM([1.0 2.0]), 2, 2)) == [1.0 2.0 1.0 2.0; 1.0 2.0 1.0 2.0]
    @test Vector{Float64}(C.linspace(DM(0.0), DM(1.0), 5)) == [0.0, 0.25, 0.5, 0.75, 1.0]
end

@testset "blockcat / veccat" begin
    # ergonomic blockcat (matrix-of-blocks) restores DM type.
    B = blockcat([[DM(1.0), DM(2.0)], [DM(3.0), DM(4.0)]])
    @test Matrix{Float64}(B) == [1.0 2.0; 3.0 4.0]
    # veccat stacks the column-vectorizations end-to-end (returns MX -> evalf).
    @test Vector{Float64}(C.evalf(C.veccat([DM([1.0, 2.0]), DM([3.0])]))) == [1.0, 2.0, 3.0]
end

# ---- type predicates: is_equal / is_op ----
@testset "is_equal / is_op predicates" begin
    x = sym(SX, "x")
    @test C.is_equal(x, x)                       # same symbol
    @test C.is_equal(C.sin(x), C.sin(x), 10)     # structural equality up to depth 10
    @test C.is_equal(C.sin(x), C.cos(x), 10) == false
    @test C.is_op(C.sin(x), C.OP_SIN)            # node opcode introspection
    @test C.is_op(C.sin(x), C.OP_COS) == false
    y = sym(SX, "y")
    @test C.is_op(x + y, C.OP_ADD)               # x+x folds to 2*x, so use x+y
end

# ---- string / printing helpers ----
@testset "str / printing" begin
    x = sym(SX, "x")
    e = C.sin(x) + x * x
    s = C.str(e)
    @test s isa AbstractString && !isempty(s)
    @test occursin("sin", s)
    # str of a Sparsity doesn't crash (regression #179 territory).
    sp = C.dense(C.Sparsity, 3, 3)
    @test C.str(sp) isa AbstractString
    # str of a Function.
    f = C.Function("f", [x], [e])
    @test occursin("f", C.str(f))
end

# ---- function call schemes (mirror misc.js test_getscheme) ----
@testset "named IO call scheme" begin
    x = sym(SX, "x"); p = sym(SX, "p")
    # 5-arg named-IO ctor needs the trailing opts dict (no defaulted overload).
    F = C.Function("F", [x, p], [x + p, x * x], ["x", "p"], ["f", "g"], Dict{String,Any}())
    fc = F(Dict("x" => DM(3.0), "p" => DM(4.0)))
    @test Float64(fc["f"]) == 7.0
    @test Float64(fc["g"]) == 9.0
    @test C.name_in(F) == ["x", "p"]
    @test C.name_out(F) == ["f", "g"]
end

# ---- DM copy independence (mirror misc.js test_copyconstr_norefcount) ----
@testset "DM copy is independent" begin
    x = C.ones(C.DM, 2, 3)
    y = C.DM(x)
    @test Vector{Float64}(x) == Vector{Float64}(y)
    @test length(C.nonzeros(x)) == 6
    @test all(Vector{Float64}(y) .== 1.0)
end

# ---- invalid-Sparsity DM construction throws, doesn't segfault (test_sanity) ----
@testset "invalid Sparsity DM throws" begin
    ok = DM(C.Sparsity(4, 3, [0, 2, 2, 3], [1, 2, 1]), DM([0.738, 0.39, 0.99]))
    @test C.nnz(ok) == 3
    @test_throws Exception DM(C.Sparsity(4, 4, [0, 2, 2, 3], [1, 2, 1]), DM([0.738, 0.39, 0.99]))
end

# ---- linear-algebra value functions with analytic checks ----
@testset "dot / norms / trace / det / inv" begin
    a = DM([1.0, 2.0, 3.0]); b = DM([4.0, 5.0, 6.0])
    @test Float64(C.dot(a, b)) == 32.0                       # 4+10+18
    @test Float64(C.norm_1(a)) == 6.0                        # |1|+|2|+|3|
    @test Float64(C.norm_2(DM([3.0, 4.0]))) == 5.0           # sqrt(9+16)
    @test Float64(C.norm_inf(DM([-7.0, 3.0, 5.0]))) == 7.0   # max |.|
    @test abs(Float64(C.norm_fro(DM([1.0 2.0; 2.0 3.0]))) - sqrt(18.0)) < 1e-12

    A = DM([2.0 1.0; 1.0 3.0])
    @test Float64(C.trace(A)) == 5.0                         # 2+3
    @test Float64(C.det(A)) == 5.0                           # 2*3 - 1*1
    @test Matrix{Float64}(C.inv(A)) == [0.6 -0.2; -0.2 0.4]  # 1/5 * [3 -1; -1 2]
    # inv(A) * A == I.
    @test maximum(abs.(Matrix{Float64}(C.mtimes(C.inv(A), A)) .- [1.0 0.0; 0.0 1.0])) < 1e-12
    # det of a 3x3 with a known value.
    @test abs(Float64(C.det(DM([1.0 2.0 3.0; 0.0 1.0 4.0; 5.0 6.0 0.0]))) - 1.0) < 1e-10
end

# ---- simplify doesn't crash and preserves value ----
@testset "simplify preserves value" begin
    x = sym(SX, "x")
    s = C.simplify(x * x + 0.0)
    @test s isa SX
    f = C.Function("f", [x], [s])
    @test abs(Float64(f(DM(3.0))) - 9.0) < 1e-12             # x^2 at x=3
end
