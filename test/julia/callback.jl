# callback.jl -- Julia-defined Callbacks (eval, derivative hooks) and opti.callback.

const CN = CasADiNative

struct CbSquare <: CasadiCallback end
CN.cb_eval(::CbSquare, a) = [a[1]^2]

# Exact Jacobian as an MX Function
struct CbSquareJac <: CasadiCallback end
CN.cb_eval(::CbSquareJac, a) = [a[1]^2]
CN.cb_has_jacobian(::CbSquareJac) = true
function CN.cb_get_jacobian(::CbSquareJac, name, inames, onames, opts)
    u = sym(MX, "u"); o = sym(MX, "o")
    C.Function(name, [u, o], [2 * u], inames, onames, opts)
end

# Exact Jacobian as another Callback (kept alive by the bindings)
struct CbTwice <: CasadiCallback end
CN.cb_n_in(::CbTwice) = 2
CN.cb_eval(::CbTwice, a) = [2 * a[1]]
struct CbSquareCbJac <: CasadiCallback end
CN.cb_eval(::CbSquareCbJac, a) = [a[1]^2]
CN.cb_has_jacobian(::CbSquareCbJac) = true
CN.cb_get_jacobian(::CbSquareCbJac, name, inames, onames, opts) = callback(CbTwice(); name = name)

# Forward and reverse sensitivities for one direction
struct CbSquareFwdRev <: CasadiCallback end
CN.cb_eval(::CbSquareFwdRev, a) = [a[1]^2]
CN.cb_has_forward(::CbSquareFwdRev, nfwd) = nfwd == 1
CN.cb_has_reverse(::CbSquareFwdRev, nadj) = nadj == 1
function _cb_sens(name, inames, onames, opts)
    u = sym(MX, "u"); o = sym(MX, "o"); s = sym(MX, "s")
    C.Function(name, [u, o, s], [2 * u * s], inames, onames, opts)
end
CN.cb_get_forward(::CbSquareFwdRev, nfwd, name, inames, onames, opts) = _cb_sens(name, inames, onames, opts)
CN.cb_get_reverse(::CbSquareFwdRev, nadj, name, inames, onames, opts) = _cb_sens(name, inames, onames, opts)

# Vector-valued, with a declared Jacobian sparsity
struct CbDiag <: CasadiCallback end
CN.cb_sparsity_in(::CbDiag, i) = C.dense(Sparsity, 3, 1)
CN.cb_sparsity_out(::CbDiag, i) = C.dense(Sparsity, 3, 1)
CN.cb_eval(::CbDiag, a) = [a[1] .^ 2]
CN.cb_has_jac_sparsity(::CbDiag, oind, iind) = true
CN.cb_get_jac_sparsity(::CbDiag, oind, iind, symmetric) = C.diag(Sparsity, 3)

mutable struct CbCounted <: CasadiCallback
    n::Int
end
CN.cb_eval(::CbCounted, a) = [a[1]]
CN.cb_init(c::CbCounted) = (c.n += 1; nothing)
CN.cb_finalize(c::CbCounted) = (c.n += 10; nothing)

# nlpsol iteration_callback: sees (x, f, g, lam_x, lam_g, lam_p) every iteration
mutable struct CbIter <: CasadiCallback
    n::Int
end
CN.cb_n_in(::CbIter) = 6
CN.cb_name_in(::CbIter, i) = C.nlpsol_out(i - 1)
CN.cb_sparsity_in(::CbIter, i) = i in (1, 4) ? C.dense(Sparsity, 2, 1) : i == 2 ? C.dense(Sparsity, 1, 1) : Sparsity(0, 0)
CN.cb_eval(c::CbIter, a) = (c.n += 1; [DM(0)])

const _ipopt_quiet = Dict("print_time" => false, "ipopt" => Dict{String,Any}("print_level" => 0))

@testset "Callback: numeric and symbolic calls" begin
    cb = callback(CbSquare(); name = "sq")
    @test Float64(cb(3)) == 9.0
    z = sym(MX, "z")
    @test cb(z) isa MX
    F = C.Function("F", [z], [cb(z) + 1])
    @test Float64(F(3)) == 10.0
end

@testset "Callback: derivative hooks" begin
    J = C.Function(callback(CbSquareJac(); name = "sqj")).jacobian()
    @test Float64(J(3, 9)) == 6.0
    cbc = callback(CbSquareCbJac(); name = "sqc")
    GC.gc(); GC.gc()
    @test Float64(C.Function(cbc).jacobian()(4, 16)) == 8.0
    G = C.Function(callback(CbSquareFwdRev(); name = "sqfr"))
    @test Float64(G.forward(1)(3, 9, 1)) == 6.0
    @test Float64(G.reverse(1)(3, 9, 1)) == 6.0
    D = C.Function(callback(CbDiag(); name = "diag", enable_fd = true))
    @test C.is_equal(D.jac_sparsity(0, 0), C.diag(Sparsity, 3))
    c = CbCounted(0)
    callback(c; name = "counted")
    @test c.n == 11
end

HAS_IPOPT && @testset "Callback inside nlpsol" begin
    z = sym(MX, "z")
    for cb in (callback(CbSquare(); name = "sq_fd", enable_fd = true),
               callback(CbSquareJac(); name = "sq_jac"))
        S = nlpsol("S", "ipopt", Dict("x" => z, "f" => (cb(z) - 2)^2),
                   merge(_ipopt_quiet, Dict("ipopt" => Dict{String,Any}("print_level" => 0,
                       "hessian_approximation" => "limited-memory"))))
        @test Base.abs(Float64(S(x0 = 1)["x"]) - sqrt(2)) < 1e-6
    end
end

HAS_IPOPT && @testset "nlpsol iteration_callback" begin
    x = sym(SX, "x", 2)
    c = CbIter(0)
    S = nlpsol("S", "ipopt", Dict("x" => x, "f" => (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2),
               merge(_ipopt_quiet, Dict("iteration_callback" => callback(c; name = "iter"))))
    S(x0 = [0, 0])
    @test c.n == Int(S.stats()["iter_count"]) + 1
end

HAS_IPOPT && @testset "opti.callback" begin
    opti = Opti()
    x = opti.variable(); y = opti.variable()
    opti.minimize((1 - x)^2 + 100 * (y - x^2)^2)
    opti.solver("ipopt", Dict("print_time" => false), Dict("print_level" => 0))
    its = Int[]; xs = Float64[]
    opti.callback(i -> (push!(its, i); push!(xs, opti.debug.value(x))))
    opti.solve()
    @test its[1:3] == [0, 1, 2]
    @test Base.abs(xs[end] - 1) < 1e-6
    opti.callback()
    empty!(its)
    opti.solve()
    @test isempty(its)
    callback(opti, i -> push!(its, i))
    opti.solve()
    @test !isempty(its)
    opti.callback(i -> error("boom"))
    @test_throws Exception opti.solve()
end
