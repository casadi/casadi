# SX symbolics, operators, indexing, derivatives, Function construction.

@testset "SX symbolics" begin
    x = sym(SX, "x")
    @test x isa SX
    @test Base.string(x) == "x"
    y = sym(SX, "y", 2, 2)
    @test Base.size(y) == (2, 2)
    @test C.is_symbolic(y)
    z = sin(x + x)
    @test Base.string(z) == "sin((2.*x))"
    @test Base.string(C.substitute(z, x, sym(SX, "q"))) == "sin((2.*q))"
end

@testset "operators and promotion" begin
    x = sym(SX, "x")
    @test Base.string(x + 1) == "(x+1)"
    @test Base.string(2 * x) == "(2.*x)"
    @test Base.string((1 - x)^2) == "sq((1-x))"
    d = DM(2.0) + DM(3.0)
    @test Float64(d) == 5.0
    @test Float64(DM(4.0) * DM(5.0)) == 20.0
    @test Float64(DM(7.0) + 1) == 8.0
end

@testset "mtimes vs broadcast" begin
    A = sym(SX, "A", 2, 2)
    b = sym(SX, "b", 2, 1)
    @test Base.size(A * b) == (2, 1)        # * is matrix product
    @test Base.size(A .* A) == (2, 2)       # .* is elementwise
    @test Base.occursin("sq", Base.string(A .* A))
end

@testset "indexing (1-based)" begin
    v = DM([10.0, 20.0, 30.0, 40.0, 50.0])
    @test Float64(v[2]) == 20.0
    @test Base.string(v[2:4]) == "[20, 30, 40]"
    @test Base.string(v[1:2:5]) == "[10, 30, 50]"
    @test Float64(v[end]) == 50.0
    m = sym(SX, "m", 2, 3)
    @test Base.size(m[1, :]) == (1, 3)
    @test Base.size(m[:, 2]) == (2, 1)
    w = DM([0.0, 0.0, 0.0])
    w[2] = 9.0
    @test Base.string(w) == "[0, 9, 0]"
end

@testset "derivatives" begin
    x = sym(SX, "x")
    @test Base.string(jacobian(sin(2 * x), x)) == "(2.*cos((2.*x)))"
    @test Base.string(jacobian(x^2, x)) == "(2.*x)"
end

@testset "Function construction and evaluation" begin
    x = sym(SX, "x")
    f = C.Function("f", [x], [sin(x + x)], Dict{String,Any}())
    @test C.n_in(f) == 1
    @test C.n_out(f) == 1
    r = f(0.5)
    @test r isa DM
    @test Base.abs(Float64(r) - Base.sin(1.0)) < 1e-12
end

@testset "exceptions" begin
    x = sym(SX, "x")
    # jacobian of incompatible shapes raises a casadi error, not a crash
    @test_throws C.SwigError C.Function("f", [x], [sym(SX, "y")], Dict{String,Any}())
end

@testset "vectors and dicts across the boundary" begin
    x = sym(SX, "x"); y = sym(SX, "y")
    w = vertcat(x, y)
    @test Base.size(w) == (2, 1)
    w2 = vertcat([x, y])
    @test Base.size(w2) == (2, 1)
    f = C.Function("f", [x], [x + 1], Dict{String,Any}())
    @test C.name_out(f, 0) == "o0"
end

@testset "comparison constraints" begin
    x = sym(SX, "x")
    @test Base.string(x <= 1) == "(x<=1)"
    @test Base.string(x == 0) == "(x==0)"
end

@testset "GC stress" begin
    for i in 1:2000
        t = sym(SX, "t")
        sin(t + DM(Float64(i)))
    end
    GC.gc()
    @test true
end
