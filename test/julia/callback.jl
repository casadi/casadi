# callback.jl -- Julia-defined Callbacks (eval, derivative hooks) and opti.callback.

const CN = CasADiNative

struct CbSquare <: Callback end
CN.eval_(::CbSquare, a) = [a[1]^2]

# Exact Jacobian as an MX Function
struct CbSquareJac <: Callback end
CN.eval_(::CbSquareJac, a) = [a[1]^2]
CN.has_jacobian(::CbSquareJac) = true
function CN.get_jacobian(::CbSquareJac, name, inames, onames, opts)
    u = sym(MX, "u"); o = sym(MX, "o")
    C.Function(name, [u, o], [2 * u], inames, onames, opts)
end

# Exact Jacobian as another Callback, kept alive by its owner (as in Python)
struct CbTwice <: Callback end
CN.get_n_in(::CbTwice) = 2
CN.eval_(::CbTwice, a) = [2 * a[1]]
mutable struct CbSquareCbJac <: Callback
    jac::Any
end
CN.eval_(::CbSquareCbJac, a) = [a[1]^2]
CN.has_jacobian(::CbSquareCbJac) = true
CN.get_jacobian(c::CbSquareCbJac, name, inames, onames, opts) = (c.jac = callback(CbTwice(), name))

# Forward and reverse sensitivities for one direction
struct CbSquareFwdRev <: Callback end
CN.eval_(::CbSquareFwdRev, a) = [a[1]^2]
CN.has_forward(::CbSquareFwdRev, nfwd) = nfwd == 1
CN.has_reverse(::CbSquareFwdRev, nadj) = nadj == 1
function _cb_sens(name, inames, onames, opts)
    u = sym(MX, "u"); o = sym(MX, "o"); s = sym(MX, "s")
    C.Function(name, [u, o, s], [2 * u * s], inames, onames, opts)
end
CN.get_forward(::CbSquareFwdRev, nfwd, name, inames, onames, opts) = _cb_sens(name, inames, onames, opts)
CN.get_reverse(::CbSquareFwdRev, nadj, name, inames, onames, opts) = _cb_sens(name, inames, onames, opts)

# Vector-valued, with a declared Jacobian sparsity
struct CbDiag <: Callback end
CN.get_sparsity_in(::CbDiag, i) = C.dense(Sparsity, 3, 1)
CN.get_sparsity_out(::CbDiag, i) = C.dense(Sparsity, 3, 1)
CN.eval_(::CbDiag, a) = [a[1] .^ 2]
CN.has_jac_sparsity(::CbDiag, oind, iind) = true
CN.get_jac_sparsity(::CbDiag, oind, iind, symmetric) = C.diag(Sparsity, 3)

mutable struct CbCounted <: Callback
    n::Int
end
CN.eval_(::CbCounted, a) = [a[1]]
CN.init(c::CbCounted) = (c.n += 1; nothing)
CN.finalize(c::CbCounted) = (c.n += 10; nothing)

# nlpsol iteration_callback: sees (x, f, g, lam_x, lam_g, lam_p) every iteration
mutable struct CbIter <: Callback
    n::Int
end
CN.get_n_in(::CbIter) = 6
CN.get_name_in(::CbIter, i) = C.nlpsol_out(i)
CN.get_sparsity_in(::CbIter, i) = i in (0, 3) ? C.dense(Sparsity, 2, 1) : i == 1 ? C.dense(Sparsity, 1, 1) : Sparsity(0, 0)
CN.eval_(c::CbIter, a) = (c.n += 1; [DM(0)])

const _ipopt_quiet = Dict("print_time" => false, "ipopt" => Dict{String,Any}("print_level" => 0))

@testset "Callback: numeric and symbolic calls" begin
    cb = callback(CbSquare(), "sq")
    @test Float64(cb(3)) == 9.0
    z = sym(MX, "z")
    @test cb(z) isa MX
    F = C.Function("F", [z], [cb(z) + 1])
    @test Float64(F(3)) == 10.0
end

@testset "Callback: base-class methods" begin
    cb = callback(CbSquare(), "sqs")
    @test occursin("sqs", C.str(cb))
    @test !C.is_null(cb)
end

struct CbNamed <: Callback end
CN.get_n_in(::CbNamed) = 2
CN.get_sparsity_in(::CbNamed, i) = i == 0 ? C.dense(Sparsity, 2, 1) : C.dense(Sparsity, 1, 1)
CN.eval_(::CbNamed, a) = [sum(a[1]) * a[2]]

@testset "Callback: 0-based hooks, i0/o0 default names" begin
    cb = callback(CbNamed(), "named")
    f = C.Function(cb)
    @test C.name_in(f) == ["i0", "i1"] && C.name_out(f) == ["o0"]
    @test C.size_in(f, 0) == (2, 1) && GC.@preserve(cb, Float64(f([1, 2], 3))) == 9
    @test C.name_in(C.Function(callback(CbSquare(), "sqn"))) == ["i0"]
end

# uses_output and buffer evaluation
struct CbBuffer <: Callback end
CN.get_n_in(::CbBuffer) = 2
CN.get_sparsity_in(::CbBuffer, i) = C.dense(Sparsity, 2, 1)
CN.get_sparsity_out(::CbBuffer, i) = C.dense(Sparsity, 2, 1)
CN.has_eval_buffer(::CbBuffer) = true
CN.eval_buffer(::CbBuffer, arg, res) = (res[1] .= arg[1] .* arg[2]; 0)
CN.uses_output(::CbBuffer) = true

@testset "Callback: positional options, eval_buffer, uses_output, no copies" begin
    cb = callback(CbBuffer(), "buf", Dict("enable_fd" => true))
    f = C.Function(cb)
    @test GC.@preserve(cb, Vector(f([1, 2], [3, 4]))) == [3, 8]
    @test C.uses_output(f)
    @test_throws MethodError C.CallbackFunction(callback(CbSquare(), "nocopy"))
end

# hook returns of any Julia number type convert losslessly; others raise
struct CbTyped <: Callback end
CN.get_n_in(::CbTyped) = Int32(2)
CN.get_n_out(::CbTyped) = true
CN.uses_output(::CbTyped) = 1
CN.get_sparsity_in(::CbTyped, i) = C.dense(Sparsity, 2, 1)
CN.get_sparsity_out(::CbTyped, i) = C.dense(Sparsity, 2, 1)
CN.has_eval_buffer(::CbTyped) = true
CN.eval_buffer(::CbTyped, arg, res) = (res[1] .= arg[1] .+ arg[2]; UInt8(0))
struct CbBadN <: Callback end
CN.get_n_in(::CbBadN) = 1.5
CN.eval_(::CbBadN, a) = [a[1]]
struct CbBadRet <: Callback end
CN.has_eval_buffer(::CbBadRet) = true
CN.eval_buffer(::CbBadRet, arg, res) = "zero"

@testset "Callback: hook return types" begin
    cb = callback(CbTyped(), "typed")
    f = C.Function(cb)
    @test (C.n_in(f), C.n_out(f)) == (2, 1)
    @test GC.@preserve(cb, Vector(f([1, 2], [3, 4]))) == [4, 6]
    @test C.uses_output(f)
    @test_throws CN.SwigError callback(CbBadN(), "badn")
    g = callback(CbBadRet(), "badret")
    @test_throws CN.SwigError g(1)
end

@testset "Callback: derivative hooks" begin
    J = C.Function(callback(CbSquareJac(), "sqj")).jacobian()
    @test Float64(J(3, 9)) == 6.0
    cbc = callback(CbSquareCbJac(nothing), "sqc")
    Jc = C.Function(cbc).jacobian()
    GC.gc(); GC.gc()
    @test Float64(Jc(4, 16)) == 8.0
    cbg = callback(CbSquareFwdRev(), "sqfr")
    G = C.Function(cbg)
    GC.@preserve cbg begin
        @test Float64(G.forward(1)(3, 9, 1)) == 6.0
        @test Float64(G.reverse(1)(3, 9, 1)) == 6.0
    end
    D = C.Function(callback(CbDiag(), "diag", Dict("enable_fd" => true)))
    @test C.is_equal(D.jac_sparsity(0, 0), C.diag(Sparsity, 3))
    c = CbCounted(0)
    callback(c, "counted")
    @test c.n == 11
end

# get_factory, and casadi's defaults for whatever is not defined
mutable struct CbFactory <: Callback
    seen::Vector{String}
end
CN.get_n_in(::CbFactory) = 2
CN.eval_(::CbFactory, a) = [a[1] * a[2]]
function CN.get_factory(c::CbFactory, name, s_in, s_out, aux, opts)
    push!(c.seen, name)
    u = sym(MX, "u"); v = sym(MX, "v")
    C.Function(name, [u, v], [u + v], s_in, s_out)
end

@testset "Callback: every virtual, casadi defaults" begin
    c = CbFactory(String[])
    cb = callback(c, "fac")
    @test C.name_in(cb) == ["i0", "i1"] && C.name_out(cb) == ["o0"]
    @test C.size_in(cb, 1) == (1, 1) && !C.has_jacobian(cb)
    F = GC.@preserve cb C.factory(cb, "sum", ["i0", "i1"], ["o0"])
    @test c.seen == ["sum"] && Float64(F(2, 3)) == 5.0
    @test CN.get_n_out(c) === C._swig_director_no_override
end

# a Function whose Callback proxy nobody keeps
@noinline _dropped_callback() = C.Function(callback(CbSquare(), "dropped"))

@testset "Callback: lifetime is the user's" begin
    f = _dropped_callback()
    GC.gc(); GC.gc(); GC.gc()
    err = try f(2); nothing catch e; e end
    @test err isa CN.SwigError && occursin("Callback object has been deleted", err.msg)
    c = CbSquareCbJac(nothing)
    J = C.Function(callback(c, "outer")).jacobian()
    c.jac = nothing
    GC.gc(); GC.gc(); GC.gc()
    err = try J(4, 16); nothing catch e; e end
    @test err isa CN.SwigError && occursin("Callback object has been deleted", err.msg)
end

HAS_IPOPT && @testset "Callback inside nlpsol" begin
    z = sym(MX, "z")
    for cb in (callback(CbSquare(), "sq_fd", Dict("enable_fd" => true)),
               callback(CbSquareJac(), "sq_jac"))
        S = nlpsol("S", "ipopt", Dict("x" => z, "f" => (cb(z) - 2)^2),
                   merge(_ipopt_quiet, Dict("ipopt" => Dict{String,Any}("print_level" => 0,
                       "hessian_approximation" => "limited-memory"))))
        @test Base.abs(GC.@preserve(cb, Float64(S(x0 = 1)["x"])) - sqrt(2)) < 1e-6
    end
end

HAS_IPOPT && @testset "nlpsol iteration_callback" begin
    x = sym(SX, "x", 2)
    c = CbIter(0)
    cb = callback(c, "iter")
    S = nlpsol("S", "ipopt", Dict("x" => x, "f" => (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2),
               merge(_ipopt_quiet, Dict("iteration_callback" => cb)))
    GC.@preserve cb S(x0 = [0, 0])
    @test c.n == Int(S.stats()["iter_count"]) + 1
end

# any callable works as an opti.callback
struct IterLog
    its::Vector{Int}
end
(l::IterLog)(i) = push!(l.its, i)

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
    log = IterLog(Int[])
    opti.callback(log)
    GC.gc(); GC.gc()
    opti.solve()
    @test log.its[1:2] == [0, 1]
    opti.callback(i -> error("boom"))
    @test_throws Exception opti.solve()
end

# casadi's worker threads call into Julia; the GC may run meanwhile
struct CbAlloc <: Callback end
CN.eval_(::CbAlloc, a) = (GC.gc(false); [a[1]^2])
@testset "Callback: evaluated by casadi worker threads" begin
    cb = callback(CbAlloc(), "sq")
    F = C.Function(cb).map(8, "thread", 4)
    GC.@preserve cb @test Base.Matrix(F(DM(Float64[1 2 3 4 5 6 7 8]))) == [1 4 9 16 25 36 49 64]
end
