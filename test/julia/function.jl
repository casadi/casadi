# function.jl -- capability port of test/javascript/function.js (the 709-line
# twin) + deeper Julia-only coverage: codegen->gcc->external round-trip,
# map/mapaccum/fold, factory, wrap, named-IO Dict & method-style calls.

const _FN_OD = Dict{String,Any}()

@testset "construct from SX/MX, call conventions, introspection" begin
    x = sym(SX, "x", 2)
    f = C.Function("f", [x], [x])
    # functor call (n_out=1 -> bare result)
    @test Vector{Float64}(f(DM([1.0, 2.0]))) == [1.0, 2.0]
    @test C.n_in(f) == 1
    @test C.n_out(f) == 1
    @test C.name(f) == "f"
    @test C.name_in(f) == ["i0"]
    @test C.name_out(f) == ["o0"]
    # C.call always returns a vector
    r = C.call(f, [DM([3.0, 4.0])])
    @test r isa AbstractVector && length(r) == 1
    @test Vector{Float64}(r[1]) == [3.0, 4.0]
end

@testset "multi-IO functor returns a tuple" begin
    x1 = sym(MX, "x", 2); y1 = sym(MX, "y")
    x2 = sym(MX, "x", 2); y2 = sym(MX, "y")
    p = C.Function("p", [x1, y1, x2, y2],
        [C.sin(x1) + y1, C.sin(x2) + y2])
    out = p(DM([4.0, 5.0]), DM(3.0), DM([5.0, 7.0]), DM(8.0))
    @test out isa Tuple && length(out) == 2
    @test maximum(abs.(Vector{Float64}(out[1]) .- (sin.([4.0, 5.0]) .+ 3.0))) < 1e-12
    @test maximum(abs.(Vector{Float64}(out[2]) .- (sin.([5.0, 7.0]) .+ 8.0))) < 1e-12
end

@testset "named IO: Dict-style call and method-style" begin
    x = sym(SX, "x")
    f = C.Function("f", [x], [C.times(x, x), x], ["i0"], ["foo", "bar"], _FN_OD)
    ret = f(Dict("i0" => DM(12.0)))
    @test Float64(ret["foo"]) == 144.0
    @test Float64(ret["bar"]) == 12.0
    # method-style introspection
    @test f.n_in() == 1
    @test f.name_out() == ["foo", "bar"]
end

@testset "sizes / sparsity introspection" begin
    x = sym(SX, "x", 3, 2)
    f = C.Function("f", [x], [x])
    @test C.size_in(f, 0) == (3, 2)
    @test C.size_out(f, 0) == (3, 2)
    @test C.size1_in(f, 0) == 3
    @test C.size2_in(f, 0) == 2
    @test C.nnz_in(f, 0) == 6
    s = C.sparsity_in(f, 0)
    @test (C.size1(s), C.size2(s)) == (3, 2)
end

@testset "sx_in / mx_in introspection (vector returns)" begin
    xs = sym(SX, "x"); ys = sym(SX, "y")
    fs = C.Function("fs", [xs, ys], [C.times(xs, ys)])
    si = C.sx_in(fs)
    @test length(si) == 2 && all(v -> v isa SX, si)
    xm = sym(MX, "x")
    fm = C.Function("fm", [xm], [C.times(xm, xm)])
    mi = C.mx_in(fm)
    @test length(mi) == 1 && mi[1] isa MX
end

@testset "expand: MX-built -> SX, numeric match" begin
    X = sym(MX, "X")
    f = C.Function("f", [X], [C.times(X, X) + 3.0])
    @test C.is_a(f, "MXFunction")
    g = C.expand(f)
    @test C.is_a(g, "SXFunction")
    for v in [0.0, 1.0, -2.5, 7.0]
        @test abs(Float64(f(DM(v))) - Float64(g(DM(v)))) < 1e-12
    end
end

@testset "wrap" begin
    x = sym(SX, "x")
    f = C.Function("f", [x], [C.times(x, x)])
    w = C.wrap(f)
    @test w isa C.Function
    @test abs(Float64(w(DM(4.0))) - 16.0) < 1e-12
end

@testset "factory: derive a gradient function" begin
    x = sym(MX, "x")
    f = C.Function("f", [x], [C.times(x, x)], ["x"], ["y"], _FN_OD)
    g = factory(f, "dfdx", ["x"], ["grad:y:x"])
    @test abs(Float64(g(DM(5.0))) - 10.0) < 1e-9      # d(x^2)/dx at 5
end

@testset "map: parallel evaluation" begin
    a = sym(SX, "a")
    g = C.Function("g", [a], [C.times(a, a)])
    m = C.map(g, 3)
    @test Vector{Float64}(m(DM([2.0 3.0 4.0]))) == [4.0, 9.0, 16.0]
end

@testset "mapaccum: running accumulation" begin
    a = sym(SX, "a")
    g = C.Function("g", [a], [a + 1.0])
    m = C.mapaccum(g, "acc", 4, _FN_OD)
    @test Vector{Float64}(m(DM(0.0))) == [1.0, 2.0, 3.0, 4.0]
end

@testset "fold: reduce to final accumulator" begin
    a = sym(SX, "a")
    g = C.Function("g", [a], [a + 1.0])
    fo = C.fold(g, 4, _FN_OD)
    @test Float64(fo(DM(0.0))) == 4.0
end

@testset "instruction-level introspection" begin
    x = sym(SX, "x")
    f1 = C.Function("f1", [x], [C.times(x, x)], Dict{String,Any}("cse" => false))
    @test C.n_instructions(f1) > 0
    @test C.is_a(f1, "SXFunction")

    # CSE must STRICTLY reduce instructions: build the same subgraph twice as
    # independent MX nodes; cse merges the duplicate sin*cos so the second copy
    # (and its sin/cos) collapses away.  A no-op `<=` would pass on equal counts.
    xm = sym(MX, "x")
    sub1 = C.times(C.sin(xm), C.cos(xm))
    sub2 = C.times(C.sin(xm), C.cos(xm))                 # structurally equal, distinct
    g1 = C.Function("g1", [xm], [sub1 + sub2], Dict{String,Any}("cse" => false))
    g2 = C.Function("g2", [xm], [sub1 + sub2], Dict{String,Any}("cse" => true))
    @test C.n_instructions(g2) < C.n_instructions(g1)    # strict reduction
    # same numeric result either way (cse changes the graph, not the value)
    @test abs(Float64(g1(DM(0.7))) - Float64(g2(DM(0.7)))) < 1e-12
end

@testset "non-symbolic input errors" begin
    x = sym(SX, "x")
    y = C.times(x, x)                                 # not a leaf symbol
    @test_throws Exception C.Function("f", [y], [x])
end

@testset "default Function is null; stats throws" begin
    f = C.Function()
    @test C.is_null(f)
    @test_throws Exception C.stats(f, true)
end

@testset "scalar broadcast on vector input" begin
    x = sym(MX, "x", 3)
    f = C.Function("f", [x], [C.times(MX(2.0), x)])
    @test Vector{Float64}(f(DM(5.0))) == [10.0, 10.0, 10.0]   # scalar 5 broadcasts
end

@testset "codegen -> gcc -> external round-trip" begin
    d = mktempdir()
    cd(d) do
        x = sym(SX, "x")
        f = C.Function("extf", [x], [2.0 * x + 1.0])
        C.generate(f, "extf.c", _FN_OD)
        @test isfile("extf.c")
        run(`gcc -fPIC -shared -O1 extf.c -o extf.so`)
        ext = C.external("extf", "./extf.so", _FN_OD)
        @test ext isa C.Function
        @test abs(Float64(ext(DM(3.0))) - 7.0) < 1e-12        # 2*3+1
        @test abs(Float64(ext(DM(-2.0))) - (-3.0)) < 1e-12
        # external function matches the original numerically
        for v in [0.0, 5.0, 10.5]
            @test abs(Float64(ext(DM(v))) - Float64(f(DM(v)))) < 1e-12
        end
    end
end

@testset "CodeGenerator: multi-function dump" begin
    x = sym(SX, "x")
    f = C.Function("f", [x], [C.sin(x)])
    g = C.Function("g", [x], [C.cos(x)])
    cg = codegen("gen")
    C.add(cg, f); C.add(cg, g)
    code = C.dump(cg)
    @test occursin("sin", code)
    @test occursin("cos", code)
end

@testset "SX/MX equivalence through a wrapped Function" begin
    builder(x) = C.sum1(C.sin(C.times(x, x)))
    xs = sym(SX, "x", 2); xm = sym(MX, "x", 2)
    fs = C.Function("fs", [xs], [builder(xs)])
    fm = C.Function("fm", [xm], [builder(xm)])
    inp = DM([1.0, 7.0])
    @test abs(Float64(fs(inp)) - Float64(fm(inp))) < 1e-10
end
