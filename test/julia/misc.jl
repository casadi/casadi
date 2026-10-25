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
    # source-tarball builds (no .git) report an empty revision; if present it must
    # be a hex sha (>= 7 chars; a full sha is 40).
    @test isempty(rev) || (occursin(r"^[0-9a-fA-F]+$", rev) && length(rev) >= 7)
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
    @test Matrix(x) == Matrix(y)
    @test length(C.nonzeros(x)) == 6
    @test all(Matrix(y) .== 1.0)
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
    # inv of an MX is the matrix inverse too
    M = sym(MX, "M", 2, 2)
    @test Matrix(C.evalf(C.substitute(inv(M), M, MX(A)))) ≈ [0.6 -0.2; -0.2 0.4]
    @test Matrix(C.evalf(C.substitute(C.inv(M), M, MX(A)))) ≈ [0.6 -0.2; -0.2 0.4]
end

# ---- simplify doesn't crash and preserves value ----
@testset "simplify preserves value" begin
    x = sym(SX, "x")
    s = C.simplify(x * x + 0.0)
    @test s isa SX
    f = C.Function("f", [x], [s])
    @test abs(Float64(f(DM(3.0))) - 9.0) < 1e-12             # x^2 at x=3
end

@testset "ca.<name> of a casadi method named like a Base function" begin
    CN = CasADiNative
    x = sym(SX, "x"); f = C.Function("f", [x], [2 * x])
    @test C.n_out(CN.map(f, 3)) == 1 && C.size2_in(CN.map(f, 3), 0) == 3
    @test C.n_in(CN.reverse(f, 1)) == 3
    @test CN.map(y -> y + 1, [1, 2]) == [2, 3]
    @test CN.max(1, 2) == 2
    @test C.nnz(CN.intersect(C.diag(Sparsity, 2), C.dense(Sparsity, 2, 2))) == 2
    @test Float64(CN.det(DM([1 2; 3 4]), "qr")) ≈ -2
    A = DM([3 1; 1 2])
    @test Float64(CN.det(C.Linsol("L", "qr", C.sparsity(A)), A)) ≈ 5
    opti = Opti(); y = opti.variable()
    for o2 in (CN.copy(opti), copy(opti))
        o2.subject_to(y >= 2)
        @test C.ng(o2) == 1 && C.ng(opti) == 0
    end
end

@testset "DaeBuilder der creates the derivative variable" begin
    dae = C.DaeBuilder("m")
    x = dae.add("x")
    @test C.str(dae.der(x)) == "der(x)"
    @test C.str(C.der(dae, x)) == "der(x)"
    @test dae.der() == ["der(x)"]
end

@testset "builders keep the type of their pieces" begin
    @test C.mtimes([DM([1 2]), DM([1 0; 0 2]), DM([1, 1])]) isa DM
    @test C.veccat([DM([1 2; 3 4]), DM(5)]) isa DM
    @test Vector(veccat(DM([1 2; 3 4]), 5)) == [1, 3, 2, 4, 5]
    @test diagcat() isa DM && size(diagcat()) == (0, 0)
    @test veccat(DM[]) isa DM && size(veccat(DM[])) == (0, 1)
    @test veccat(SX[]) isa SX
    @test Matrix(blockcat(1, 2, 3, 4)) == [1 2; 3 4]
    @test C._horzcat([DM(1), DM(2)]) isa DM
end

@testset "nested-vector, dict, hash and data-member returns are wrapped" begin
    @test C.blocksplit(DM([1 2; 3 4]), 1, 1) isa Vector{Vector{DM}}
    x = sym(SX, "x")
    @test C.str(C.forward([x^2], [x], [[SX(1)]])[1][1]) == "(2.*x)"
    @test length(C.reverse([x^2], [x], [[SX(1)]])) == 1
    @test size(C.sym(MX, "x", 1, 1, 2, 2)) == (2,) && C.sym(SX, "y", 2, 3, 2, 2)[1][1] isa SX
    X = sym(MX, "X", 3, 3)
    @test C.mapping(X[1, 2]) == reshape([3], 1, 1)
    @test C.hash(C.dense(Sparsity, 2)) == C.hash(C.dense(Sparsity, 2))
    @test CasADiNative.hash(C.dense(Sparsity, 2)) == C.hash(C.dense(Sparsity, 2))
    opti = Opti(); z = opti.variable()
    c = C.canon_expr(opti.debug, 0 <= (z <= 1))
    @test c.type == C.OPTI_DOUBLE_INEQUALITY && C.str(c.canon) == C.str(z)
end

@testset "DaeBuilder get returns a native value" begin
    dae = C.DaeBuilder("m")
    dae.add("x", Dict("start" => 1))
    @test dae.get("x") == 1.0
    dae.add("v", Dict("dimension" => [2], "start" => 3))
    @test dae.get(["x", "v"]) == [1.0, [3.0, 3.0]]
end

@testset "display is casadi's repr, print its str, as Python" begin
    x = sym(MX, "x")
    @test repr("text/plain", x + 1) == "MX((1+x))"
    @test repr("text/plain", x) == "MX(x)"
    @test !occursin("casadi", repr("text/plain", sym(MX, "y", 2, 3)'))
    @test repr("text/plain", DM([1.17708, 1e-7])) == "DM([1.17708, 1e-07])"
    @test repr("text/plain", DM([3 1.5])) == "DM([[3, 1.5]])"
    @test repr("text/plain", C.inf(DM)) == "DM(inf)"
    @test repr("text/plain", DM(3.5)) == "DM(3.5)"
    @test repr("text/plain", DM(2, 0)) == "DM(2x0)"
    @test repr("text/plain", sum(DM(3, 0))) == "DM(00)"
    @test repr("text/plain", DM(C.diag(Sparsity, 2), [1, 2])) == "DM(\n[[1, 00], \n [00, 2]])"
    @test repr("text/plain", sym(SX, "x", C.dense(Sparsity, 2, 2))) == "SX(\n[[x_0, x_2], \n [x_1, x_3]])"
    @test repr(sym(SX, "x", C.dense(Sparsity, 2, 2))) == "\n[[x_0, x_2], \n [x_1, x_3]]"
    @test repr(DM([1 2; 3 4])) == "\n[[1, 2], \n [3, 4]]"
    @test "$(DM([1, 2]))" == "[1, 2]"
end

@testset "Function display is casadi's summary" begin
    x = sym(MX, "x"); f = C.Function("f", [x], [x])
    @test repr("text/plain", f) == "Function(f:(i0)->(o0) MXFunction)"
    @test repr("text/plain", C.expand(f)) == "Function(f:(i0)->(o0) SXFunction)"
end

@testset "display of the other casadi objects" begin
    @test repr("text/plain", C.Linsol("L", "ldl", C.dense(Sparsity, 2, 2))) == "Linsol(Linear solver LinsolLdl)"
    @test startswith(repr("text/plain", C.DaeBuilder("m")), "DaeBuilder(nt = 0")
    @test startswith(repr("text/plain", Opti().debug), "Opti(Opti {")
    @test repr("text/plain", C.Slice(0, 3)) == "Slice(:3)"
    @test !occursin("Ptr", repr(C.Linsol("L", "ldl", C.dense(Sparsity, 2, 2))))
end

@testset "bound methods display their name" begin
    @test repr("text/plain", C.Linsol.load_plugin) == "Linsol.load_plugin"
    @test repr(sym(SX, "x").nnz) == "nnz bound to a SX"
    @test Opti().variable isa Base.Function
end

@testset "errors speak Julia" begin
    @test_throws DimensionMismatch [C.ones(DM, 2, 2) C.ones(DM, 3, 1)]
    @test_throws DimensionMismatch [C.ones(DM, 2, 2); C.ones(DM, 1, 3)]
    @test_throws DimensionMismatch [sym(SX, "a", 2) sym(SX, "b", 3)]
    @test size([DM(1) DM(Sparsity(0, 0)) DM(2)]) == (1, 2)
    e = try C.Function("f", [sym(SX, "x") + 1], [SX(0)]) catch err; err end
    @test !occursin("Main.", repr(e))
    @test CasADiNative.zeros(DM, C.diag(Sparsity, 2)) isa DM
    @test repr(Vector{MX}) == "Vector{MX}" && repr(MX[]) == "MX[]"
    @test repr((DM[DM(1), DM(2)], 0)) == "(DM[1, 2], 0)"
end

@testset "casadi paths set at load, as in Python" begin
    @test isdir(C.GlobalOptions.getCasadiPath())
    @test endswith(C.GlobalOptions.getCasadiIncludePath(), "include")
end

@testset "printout goes through Julia's stdout, in order" begin
    x = sym(MX, "x")
    f = C.Function("f", [x], [C.monitor(2x, "mon")])
    path, io = mktemp()
    redirect_stdout(io) do
        println("before"); f(1); println("after")
    end
    close(io)
    out = read(path, String)
    @test findfirst("before", out)[1] < findfirst("mon:", out)[1] < findfirst("after", out)[1]
end

@testset "propertynames lists methods, errors print as SwigError" begin
    x = sym(MX, "x")
    @test :jacobian in propertynames(C.Function("f", [x], [x]))
    @test :sparsity in propertynames(x) && :nz in propertynames(x)
    @test !any(n -> startswith(string(n), "_"), propertynames(C.Opti()))
    err = try C.Function("f", [x], [sym(MX, "y")]) catch e; e end
    @test err isa C.SwigError && string(typeof(err)) == "SwigError"
end
