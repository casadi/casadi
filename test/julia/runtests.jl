# test/julia/runtests.jl -- test suite for the SWIG-generated Julia bindings.
# Usage: julia runtests.jl   (needs CASADI_JL pointing at the dir holding
# CasADi.jl/casadi.jl/libcasadi_wrap.so, and casadi plugins on the search path)
using Test

const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
using .CasADiNative
const C = CasADiNative.C

const HAS_IPOPT = try
    CasADiNative.C.has_nlpsol("ipopt")
catch
    false
end
HAS_IPOPT || @info "ipopt plugin not found: solver testsets will be skipped"

@testset "casadi julia bindings" begin

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
    f = CasadiFunction("f", [x], [sin(x + x)], Dict{String,Any}())
    @test C.n_in(f) == 1
    @test C.n_out(f) == 1
    r = f(0.5)
    @test r isa DM
    @test Base.abs(Float64(r) - Base.sin(1.0)) < 1e-12
end

@testset "exceptions" begin
    x = sym(SX, "x")
    # jacobian of incompatible shapes raises a casadi error, not a crash
    @test_throws C.SwigError CasadiFunction("f", [x], [sym(SX, "y")], Dict{String,Any}())
end

@testset "vectors and dicts across the boundary" begin
    x = sym(SX, "x"); y = sym(SX, "y")
    w = vertcat(x, y)
    @test Base.size(w) == (2, 1)
    w2 = vertcat([x, y])
    @test Base.size(w2) == (2, 1)
    f = CasadiFunction("f", [x], [x + 1], Dict{String,Any}())
    @test C.name_out(f, 0) == "o0"
end

HAS_IPOPT && @testset "nlpsol rosenbrock (ipopt)" begin
    x = sym(SX, "x"); y = sym(SX, "y"); z = sym(SX, "z")
    S = nlpsol("S", "ipopt",
               Dict("x" => vertcat(x, y, z), "f" => x^2 + 100 * z^2,
                    "g" => z + (1 - x)^2 - y);
               ipopt = Dict{String,Any}("print_level" => 0), print_time = false)
    sol = S(Dict("x0" => DM([2.5, 3.0, 0.75]), "lbg" => DM(0.0), "ubg" => DM(0.0)))
    xs = sol["x"]
    @test Base.abs(Float64(xs[2]) - 1.0) < 1e-8
    @test Base.abs(Float64(xs[1])) < 1e-8
    @test Float64(sol["f"]) < 1e-10
end

HAS_IPOPT && @testset "kwarg solver call" begin
    x = sym(SX, "x"); y = sym(SX, "y"); z = sym(SX, "z")
    S = nlpsol("S", "ipopt",
               Dict("x" => vertcat(x, y, z), "f" => x^2 + 100 * z^2,
                    "g" => z + (1 - x)^2 - y);
               ipopt = Dict{String,Any}("print_level" => 0), print_time = false)
    sol = S(x0 = [2.5, 3.0, 0.75], lbg = 0, ubg = 0)   # bare numbers/arrays
    @test Base.abs(Float64(sol["x"][2]) - 1.0) < 1e-8
end

@testset "MX path" begin
    xm = sym(MX, "x", 2)
    @test xm isa MX
    fm = CasadiFunction("fm", [xm], [xm + xm], Dict{String,Any}())
    r = fm(DM([1.0, 2.0]))
    @test Vector{Float64}(r) == [2.0, 4.0]
    @test Base.size(jacobian(xm + xm, xm)) == (2, 2)
    @test Base.string(2 * sym(MX, "q")) != ""   # MX promotion
end

HAS_IPOPT && @testset "Opti stack" begin
    opti = Opti()
    x = variable(opti)
    y = variable(opti)
    minimize(opti, (x - 1)^2 + y^2)
    subject_to(opti, bounded(0, x, 0.5))   # instance-style bounded, bare numbers
    solver!(opti, "ipopt"; ipopt = Dict{String,Any}("print_level" => 0), print_time = false)
    sol = solve!(opti)
    @test Base.abs(value(sol, x) - 0.5) < 1e-6
    @test Base.abs(value(sol, y)) < 1e-6
end

@testset "comparison constraints + 2-D Opti" begin
    opti = Opti()
    X = variable(opti, 2, 3)
    @test Base.size(X[1, :]) == (1, 3)
    @test Base.size(X[:, 2]) == (2, 1)
    con = X[1, 1] == 0
    @test con isa MX
    con2 = X[2, 1] <= 1
    @test con2 isa MX
    x = sym(SX, "x")
    @test Base.string(x <= 1) == "(x<=1)"
    @test Base.string(x == 0) == "(x==0)"
end

HAS_IPOPT && @testset "double-bounded constraint (python-style chaining)" begin
    opti = Opti()
    v = variable(opti)
    minimize(opti, (v - 3)^2)
    subject_to(opti, (0 <= v) <= 1)     # (a<=b)<=c forms a double bound
    solver!(opti, "ipopt"; ipopt = Dict{String,Any}("print_level" => 0), print_time = false)
    sol = solve!(opti)
    @test Base.abs(value(sol, v) - 1.0) < 1e-6
end

@testset "GC stress" begin
    for i in 1:2000
        t = sym(SX, "t")
        sin(t + DM(Float64(i)))
    end
    GC.gc()
    @test true
end

end
