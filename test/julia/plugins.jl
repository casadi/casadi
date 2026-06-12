# plugins.jl -- capability port of test/javascript/plugins.js (+ deeper).
# Plugin enumeration / availability / loading / options-introspection across
# every plugin family, plus a trivial construct-and-solve smoke per plugin so
# "it registers" is backed by "it actually computes".  Where the JS twin only
# probes a fixed plugin list, this also exercises the meta enumeration API and
# the doc/options introspection that JS lacks.

const _PL_OD = Dict{String,Any}()

# ---- meta enumeration: plugins() / feature_list() registry strings ----
@testset "CasadiMeta plugins/features enumeration" begin
    pl = C.plugins(C.CasadiMeta)
    @test pl isa AbstractString
    @test occursin("Nlpsol::ipopt", pl)        # registered nlpsol plugin
    @test occursin("Conic::qrqp", pl)           # registered conic plugin
    @test occursin("Integrator::rk", pl)        # registered integrator plugin
    @test occursin("Rootfinder::newton", pl)    # registered rootfinder plugin
    @test occursin("Linsol::qr", pl)            # registered linsol plugin
    fl = C.feature_list(C.CasadiMeta)
    @test fl isa AbstractString
    @test occursin("interface", fl)             # human-readable feature blurbs
end

# ---- availability predicates: assert specific known-present plugins ----
@testset "has_* availability predicates" begin
    @test C.has_nlpsol("ipopt")
    @test C.has_nlpsol("sqpmethod")
    @test C.has_conic("qrqp")
    @test C.has_conic("ipqp")
    @test C.has_integrator("rk")
    @test C.has_integrator("collocation")
    @test C.has_rootfinder("newton")
    @test C.has_rootfinder("kinsol")
    @test C.has_linsol("qr")
    @test C.has_linsol("ldl")
    @test C.has_interpolant("linear")
    @test C.has_interpolant("bspline")
    # negative case: a bogus plugin name reports absent (no throw).
    @test C.has_nlpsol("nope_nope") == false
    @test C.has_conic("nope_nope") == false
    @test C.has_integrator("nope_nope") == false
    @test C.has_rootfinder("nope_nope") == false
    @test C.has_linsol("nope_nope") == false
    @test C.has_interpolant("nope_nope") == false
end

# ---- explicit plugin loading: load_* succeeds for present, throws for bogus ----
@testset "load_* plugin loading" begin
    # load is idempotent for already-registered plugins (warns "already in use").
    @test (C.load_nlpsol("ipopt"); true)
    @test (C.load_conic("qrqp"); true)
    @test (C.load_linsol("qr"); true)
    # bogus plugin: load throws cleanly (dlopen failure surfaced as an exception).
    @test_throws Exception C.load_nlpsol("nope_nope")
    @test_throws Exception C.load_conic("nope_nope")
end

# ---- options / doc introspection (JS twin has none of this) ----
@testset "options & doc introspection" begin
    # nlpsol_options(plugin) lists the option names the solver accepts.
    no = C.nlpsol_options("ipopt")
    @test no isa AbstractVector
    @test "ipopt" in no                        # ipopt's nested-options key
    # conic_options(plugin): qrqp exposes its iteration knobs.
    co = C.conic_options("qrqp")
    @test "max_iter" in co
    @test "print_iter" in co
    # nlpsol_in / nlpsol_out: the standard NLP IO scheme names.
    @test C.nlpsol_in() == ["x0", "p", "lbx", "ubx", "lbg", "ubg", "lam_x0", "lam_g0"]
    @test C.nlpsol_out() == ["x", "f", "g", "lam_x", "lam_g", "lam_p"]
    # doc_nlpsol(plugin): a non-empty documentation blob.
    doc = C.doc_nlpsol("ipopt")
    @test doc isa AbstractString && length(doc) > 100
end

# ---- bogus plugin construction throws cleanly (no segfault) ----
@testset "bogus plugin construction throws" begin
    H = DM([1.0 0.0; 0.0 1.0])
    @test_throws Exception C.conic("x", "nope_nope", Dict("h" => C.sparsity(H)), _PL_OD)
    x = sym(SX, "x")
    @test_throws Exception C.nlpsol("x", "nope_nope", Dict{String,Any}("x" => x, "f" => x * x), _PL_OD)
end

# ---- linsol plugins: A=diag([1,2,3]), x = b ./ [1,2,3] (mirrors plugins.js) ----
function _pl_linsol(plugin)
    A = DM([1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0])
    ls = C.Linsol("ls", plugin, C.sparsity(A), _PL_OD)
    C.sfact(ls, A); C.nfact(ls, A)
    x = C.solve(ls, A, DM([6.0, 8.0, 9.0]), false)
    @test maximum(abs.(Vector{Float64}(x) .- [6.0, 4.0, 3.0])) < 1e-8
end
for plugin in ["qr", "ldl", "lsqr", "symbolicqr", "csparse", "lapacklu", "lapackqr"]
    (try C.has_linsol(plugin) catch; false end) || continue
    @testset "plugins linsol $plugin" begin _pl_linsol(plugin) end
end

# tridiag needs an actual tri-banded sparsity (diag-only reads OOB) -- mirror JS.
@testset "plugins linsol tridiag" begin
    if !(try C.has_linsol("tridiag") catch; false end)
        @test_skip "tridiag linsol not loaded"
    else
        sp = C.triplet(C.Sparsity, 3, 3, [0, 1, 0, 1, 2, 1, 2], [0, 0, 1, 1, 1, 2, 2])
        A = DM(sp, DM([2.0, 1.0, 1.0, 2.0, 1.0, 1.0, 2.0]))   # tridiag [[2,1,0],[1,2,1],[0,1,2]]
        ls = C.Linsol("ls", "tridiag", sp, _PL_OD)
        C.sfact(ls, A); C.nfact(ls, A)
        x = C.solve(ls, A, DM([4.0, 8.0, 8.0]), false)         # A*[1,2,3] = [4,8,8]
        @test maximum(abs.(Vector{Float64}(x) .- [1.0, 2.0, 3.0])) < 1e-6
    end
end

# ---- integrator plugins: dx/dt = x, x(0)=1 -> x(1) = e ----
function _pl_integrator(plugin, places)
    x = sym(SX, "x")
    F = integrator("F", plugin, Dict("x" => x, "ode" => x), 0.0, 1.0)
    r = F(Dict("x0" => DM(1.0)))
    @test abs(Float64(r["xf"]) - exp(1.0)) < 10.0^(-places)
end
@testset "plugins integrator rk"          begin _pl_integrator("rk", 3) end
@testset "plugins integrator collocation" begin _pl_integrator("collocation", 4) end
for plugin in ["cvodes", "idas"]
    (try C.has_integrator(plugin) catch; false end) || continue
    @testset "plugins integrator $plugin" begin _pl_integrator(plugin, 4) end
end

# ---- interpolant plugins: linear & bspline of f(0)=0,f(1)=1,f(2)=4,f(3)=9 ----
@testset "plugins interpolant linear" begin
    f = C.interpolant("f", "linear", [[0.0, 1.0, 2.0, 3.0]], [0.0, 1.0, 4.0, 9.0], _PL_OD)
    @test abs(Float64(f(DM(1.5))) - 2.5) < 1e-9      # mid of [1,4]
    @test abs(Float64(f(DM(2.0))) - 4.0) < 1e-9      # exact at knot
end
@testset "plugins interpolant bspline" begin
    f = C.interpolant("f", "bspline", [[0.0, 1.0, 2.0, 3.0]], [0.0, 1.0, 4.0, 9.0], _PL_OD)
    @test abs(Float64(f(DM(0.0))) - 0.0) < 1e-6      # clamped left knot
    @test abs(Float64(f(DM(3.0))) - 9.0) < 1e-6      # clamped right knot
end

# ---- rootfinder plugins: x^2 - 2 = 0 from x0=1.5 -> sqrt(2) ----
function _pl_rootfinder(plugin, opts)
    x = sym(SX, "x")
    f = C.Function("f", [x], [C.times(x, x) - 2.0])
    F = C.rootfinder("F", plugin, f, opts)
    @test abs(Float64(F(DM(1.5))) - sqrt(2)) < 1e-6
end
@testset "plugins rootfinder newton"      begin _pl_rootfinder("newton", Dict{String,Any}("print_iteration" => false)) end
@testset "plugins rootfinder fast_newton" begin _pl_rootfinder("fast_newton", _PL_OD) end
@testset "plugins rootfinder kinsol" begin
    (try C.has_rootfinder("kinsol") catch; false end) ?
        _pl_rootfinder("kinsol", Dict{String,Any}("print_level" => 0)) :
        @test_skip "kinsol not loaded"
end
HAS_IPOPT && @testset "plugins rootfinder nlpsol(implicit_to_nlp)" begin
    _pl_rootfinder("nlpsol", Dict{String,Any}("nlpsol" => "ipopt",
        "nlpsol_options" => Dict{String,Any}("print_time" => false,
            "ipopt" => Dict{String,Any}("print_level" => 0, "sb" => "yes"))))
end

# ---- conic plugins: min 0.5(x^2+y^2) - 0.7x - 2.3y -> [0.7, 2.3] ----
function _pl_conic(plugin, opts)
    H = DM([1.0 0.0; 0.0 1.0]); G = DM([-0.7, -2.3])
    S = C.conic("solver", plugin, Dict("h" => C.sparsity(H)), opts)
    r = S(Dict("h" => H, "g" => G, "lbx" => DM([-Inf, -Inf]), "ubx" => DM([Inf, Inf])))
    @test maximum(abs.(Vector{Float64}(r["x"]) .- [0.7, 2.3])) < 1e-4
end
@testset "plugins conic qrqp"  begin _pl_conic("qrqp", Dict{String,Any}("print_iter" => false, "print_header" => false)) end
@testset "plugins conic ipqp"  begin _pl_conic("ipqp", Dict{String,Any}("print_iter" => false, "print_header" => false)) end
@testset "plugins conic highs" begin
    (try C.has_conic("highs") catch; false end) ?
        _pl_conic("highs", Dict{String,Any}("highs" => Dict{String,Any}("output_flag" => false))) :
        @test_skip "highs not loaded"
end
@testset "plugins conic daqp" begin
    (try C.has_conic("daqp") catch; false end) ? _pl_conic("daqp", _PL_OD) : @test_skip "daqp not loaded"
end
@testset "plugins conic qpoases" begin
    (try C.has_conic("qpoases") catch; false end) ?
        _pl_conic("qpoases", Dict{String,Any}("printLevel" => "none")) : @test_skip "qpoases not loaded"
end
HAS_IPOPT && @testset "plugins conic nlpsol(qp_to_nlp)" begin
    _pl_conic("nlpsol", Dict{String,Any}("nlpsol" => "ipopt",
        "nlpsol_options" => Dict{String,Any}("print_time" => false,
            "ipopt" => Dict{String,Any}("print_level" => 0))))
end

# ---- nlpsol plugins: f(x,y) = (x-1)^2 + (y-2)^2 -> (1, 2), f* = 0 ----
function _pl_nlpsol_sx(plugin, opts)
    xy = sym(SX, "xy", 2)
    f = (xy[1] - 1.0)^2 + (xy[2] - 2.0)^2
    S = C.nlpsol("S", plugin, Dict{String,Any}("x" => xy, "f" => f), opts)
    r = S(Dict("x0" => DM([0.0, 0.0])))
    @test maximum(abs.(Vector{Float64}(r["x"]) .- [1.0, 2.0])) < 1e-3
    @test abs(Float64(r["f"])) < 1e-4
end
@testset "plugins nlpsol sqpmethod" begin
    _pl_nlpsol_sx("sqpmethod", Dict{String,Any}("qpsol" => "qrqp",
        "qpsol_options" => Dict{String,Any}("print_iter" => false, "print_header" => false),
        "print_iteration" => false, "print_header" => false, "print_status" => false, "print_time" => false))
end
@testset "plugins nlpsol qrsqp" begin
    (try C.has_nlpsol("qrsqp") catch; false end) ?
        _pl_nlpsol_sx("qrsqp", Dict{String,Any}("print_iteration" => false, "print_header" => false, "print_time" => false)) :
        @test_skip "qrsqp not loaded"
end
@testset "plugins nlpsol feasiblesqpmethod" begin
    (try C.has_nlpsol("feasiblesqpmethod") catch; false end) ?
        _pl_nlpsol_sx("feasiblesqpmethod", Dict{String,Any}("qpsol" => "qrqp",
            "qpsol_options" => Dict{String,Any}("print_iter" => false, "print_header" => false),
            "print_iteration" => false, "print_header" => false, "print_status" => false, "print_time" => false)) :
        @test_skip "feasiblesqpmethod not loaded"
end
HAS_IPOPT && @testset "plugins nlpsol ipopt" begin
    _pl_nlpsol_sx("ipopt", Dict{String,Any}("print_time" => false,
        "ipopt" => Dict{String,Any}("print_level" => 0, "sb" => "yes")))
end

# scpgen needs an MXFunction oracle (generate_lifted is MX-only) -- mirror JS.
@testset "plugins nlpsol scpgen" begin
    if !(try C.has_nlpsol("scpgen") catch; false end)
        @test_skip "scpgen not loaded"
    else
        xy = sym(MX, "xy", 2)
        f = (xy[1] - 1.0)^2 + (xy[2] - 2.0)^2
        S = C.nlpsol("S", "scpgen", Dict{String,Any}("x" => xy, "f" => f),
            Dict{String,Any}("qpsol" => "qrqp",
                "qpsol_options" => Dict{String,Any}("print_iter" => false, "print_header" => false),
                "print_header" => false, "print_time" => false))
        r = S(Dict("x0" => DM([0.0, 0.0])))
        @test maximum(abs.(Vector{Float64}(r["x"]) .- [1.0, 2.0])) < 1e-3
        @test abs(Float64(r["f"])) < 1e-4
    end
end
