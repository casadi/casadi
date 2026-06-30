# P1-4: qpsol / conic / integrator / rootfinder / codegen convenience wrappers.

const HAS_QRQP = try C.has_conic("qrqp") catch; false end
const HAS_OSQP = try C.has_conic("osqp") catch; false end
const HAS_RK = try C.has_integrator("rk") catch; false end
const HAS_NEWTON = try C.has_rootfinder("newton") catch; false end

HAS_QRQP && @testset "qpsol (qrqp)" begin
    x = sym(SX, "x"); y = sym(SX, "y")
    f = 0.5 * (x^2 + y^2) - 2 * x - 3 * y
    Q = qpsol("Q", "qrqp", Dict("x" => vertcat(x, y), "f" => f), Dict("print_time" => false))
    sol = Q(Dict("x0" => DM([0.0, 0.0])))
    xs = Vector{Float64}(sol["x"])
    @test Base.abs(xs[1] - 2.0) < 1e-6
    @test Base.abs(xs[2] - 3.0) < 1e-6
end

HAS_OSQP && @testset "conic (osqp)" begin
    H = DM([1.0 0.0; 0.0 1.0])
    co = conic("Co", "osqp", Dict("h" => C.sparsity(H), "a" => C.Sparsity(0, 2)), Dict("print_time" => false))
    sol = co(Dict("h" => H, "g" => DM([-2.0, -3.0])))
    xs = Vector{Float64}(sol["x"])
    @test Base.abs(xs[1] - 2.0) < 1e-4
    @test Base.abs(xs[2] - 3.0) < 1e-4
end

HAS_RK && @testset "integrator (rk): xdot = -x over [0,1]" begin
    x = sym(SX, "x")
    F = integrator("F", "rk", Dict("x" => x, "ode" => -x), 0.0, 1.0)
    res = F(Dict("x0" => DM(1.0)))
    @test Base.abs(Float64(res["xf"]) - Base.exp(-1.0)) < 1e-4
end

HAS_NEWTON && @testset "rootfinder (newton): x^2 - 2 = 0" begin
    x = sym(SX, "x")
    g = C.Function("g", [x], [x^2 - 2.0], Dict{String,Any}())
    R = rootfinder("R", "newton", g)
    @test Base.abs(Float64(R(DM(1.0))) - Base.sqrt(2.0)) < 1e-9
end

@testset "codegen convenience" begin
    x = sym(SX, "x")
    f = C.Function("f", [x], [sin(x)], Dict{String,Any}())
    cg = codegen("gen")
    C.add(cg, f)
    code = C.dump(cg)
    @test Base.occursin("sin", code)
    @test cg isa CodeGenerator
end
