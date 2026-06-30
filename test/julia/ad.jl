# ad.jl -- capability port of test/javascript/ad.js (+ deeper).
# jacobian / gradient / hessian / jtimes / forward / reverse, Jacobian
# sparsity, hessian-of-Lagrangian, all checked against analytic values.

const _AD_OD = Dict{String,Any}()

@testset "jacobian: scalar SX and MX" begin
    xs = sym(SX, "x")
    fs = C.Function("fs", [xs], [jacobian(sin(xs), xs)])
    @test abs(Float64(fs(DM(1.0))) - cos(1.0)) < 1e-10
    xm = sym(MX, "x")
    fm = C.Function("fm", [xm], [jacobian(C.times(xm, xm), xm)])
    @test abs(Float64(fm(DM(3.0))) - 6.0) < 1e-10   # d(x^2)/dx = 2x
end

@testset "jacobian: vector, sparsity + values" begin
    x = sym(SX, "x", 2)
    e = vertcat(C.times(x[1], x[2]), x[1] + x[2])     # f = [x0*x1, x0+x1]
    J = jacobian(e, x)
    @test Base.size(J) == (2, 2)
    # J = [[x1, x0], [1, 1]]; at x=[3,5] column-major nonzeros [5,1,3,1]
    f = C.Function("J", [x], [J])
    @test Vector{Float64}(f(DM([3.0, 5.0]))) == [5.0, 1.0, 3.0, 1.0]
    # Diagonal Jacobian has a sparse pattern (3 nnz on a 3x3).
    y = sym(SX, "y", 3)
    Jd = jacobian(C.times(y, y), y)
    @test C.nnz(C.sparsity(Jd)) == 3
end

@testset "gradient: analytic check" begin
    x = sym(SX, "x"); y = sym(SX, "y")
    f = C.times(x, x) + C.times(C.times(SX(3.0), x), y) + C.times(y, y)
    z = vertcat(x, y)
    g = gradient(f, z)                                # [2x+3y, 3x+2y]
    @test Base.size(g) == (2, 1)
    gf = C.Function("g", [z], [g])
    @test Vector{Float64}(gf(DM([2.0, 5.0]))) == [2 * 2 + 3 * 5, 3 * 2 + 2 * 5]
end

@testset "hessian: 2nd derivatives" begin
    x = sym(SX, "x", 2)
    f = C.times(x[1], x[1]) + C.times(C.times(x[1], x[2]), x[2])  # x0^2 + x0*x1^2
    H = hessian(f, x)
    @test Base.size(H) == (2, 2)
    hf = C.Function("H", [x], [H])
    # H = [[2, 2*x1], [2*x1, 2*x0]]; at x=[3,5] -> [[2,10],[10,6]]
    @test Matrix{Float64}(hf(DM([3.0, 5.0]))) == [2.0 10.0; 10.0 6.0]
    # x^4 -> d2/dx2 = 12 x^2 ; at 2 -> 48
    s = sym(SX, "s")
    h2 = C.Function("h2", [s], [hessian(s^4, s)])
    @test abs(Float64(h2(DM(2.0))) - 48.0) < 1e-9
end

@testset "jtimes: forward jacobian-vector product" begin
    y = sym(SX, "y", 2)
    e = vertcat(C.times(y[1], y[2]), y[1] + y[2])     # J = [[y1,y0],[1,1]]
    v = sym(SX, "v", 2)
    jt = C.jtimes(e, y, v, false, _AD_OD)             # J * v
    f = C.Function("jt", [y, v], [jt])
    # at y=[3,5], v=[1,0] -> [y1*1, 1*1] = [5, 1]
    @test Vector{Float64}(f(DM([3.0, 5.0]), DM([1.0, 0.0]))) == [5.0, 1.0]
    # at v=[0,1] -> [y0, 1] = [3, 1]
    @test Vector{Float64}(f(DM([3.0, 5.0]), DM([0.0, 1.0]))) == [3.0, 1.0]
end

@testset "forward / reverse directional derivatives" begin
    x = sym(SX, "x")
    f = C.Function("f", [x], [C.times(x, x)])     # f(x) = x^2
    # forward: inputs (x, nominal_out, seed) -> 2*x*seed
    fwd = C.forward(f, 1)
    @test C.n_in(fwd) == 3
    @test abs(Float64(fwd(DM(3.0), DM(9.0), DM(1.0))) - 6.0) < 1e-10
    @test abs(Float64(fwd(DM(3.0), DM(9.0), DM(2.0))) - 12.0) < 1e-10
    # reverse: same chain rule, adjoint seed -> 2*x*seed
    rev = C.reverse(f, 1)
    @test C.n_in(rev) == 3
    @test abs(Float64(rev(DM(3.0), DM(9.0), DM(1.0))) - 6.0) < 1e-10
end

@testset "hessian of the Lagrangian via factory" begin
    x = sym(SX, "x", 2); p = sym(SX, "p")
    nlp = Dict("x" => x, "p" => p, "f" => x[1]^2 + C.times(p, x[2]^2), "g" => x[1] + x[2])
    f = C.Function("f", nlp, ["x", "p"], ["f", "g"], _AD_OD)
    H = C.factory(f, "H", ["x", "p"], ["hess:f:x:x"],
        Dict{String,Vector{String}}(), _AD_OD)
    @test C.n_out(H) == 1
    # hess(x0^2 + p*x1^2) wrt x = diag(2, 2p); at p=3 -> diag(2,6)
    @test Matrix{Float64}(H(DM([1.0, 1.0]), DM(3.0))) == [2.0 0.0; 0.0 6.0]
end

@testset "SX vs MX equivalence on shared expressions" begin
    build(x) = C.sin(x) + C.cos(x)
    xs = sym(SX, "x", 4); xm = sym(MX, "x", 4)
    fs = C.Function("fs", [xs], [build(xs)])
    fm = C.Function("fm", [xm], [build(xm)])
    inp = DM([0.0, 0.7, 1.4, 2.1])
    @test maximum(abs.(Vector{Float64}(fs(inp)) .- Vector{Float64}(fm(inp)))) < 1e-12
    # jacobians agree too
    Js = C.Function("Js", [xs], [jacobian(build(xs), xs)])
    Jm = C.Function("Jm", [xm], [jacobian(build(xm), xm)])
    @test maximum(abs.(Vector{Float64}(Js(inp)) .- Vector{Float64}(Jm(inp)))) < 1e-12
end
