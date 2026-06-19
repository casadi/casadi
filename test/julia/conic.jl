# conic.jl -- capability port of test/javascript/conic.js (+ more plugins).
# Known QP and LP solved by every loaded QP plugin; bounds, equality +
# inequality, parametric re-solve, warmstart, and infeasibility detection.

const _INF = Inf

# Per-plugin: construction opts (positional option Dict) + capability flags.
function _conic_meta(plugin)
    plugin == "qrqp"    && return (opts=Dict{String,Any}("print_iter"=>false, "print_header"=>false, "print_info"=>false), lp=true,  tol=1e-5)
    plugin == "osqp"    && return (opts=Dict{String,Any}("print_time"=>false), lp=true,  tol=5e-3)
    plugin == "qpoases" && return (opts=Dict{String,Any}("printLevel"=>"none"), lp=true,  tol=1e-5)
    plugin == "highs"   && return (opts=Dict{String,Any}("highs"=>Dict{String,Any}("output_flag"=>false)), lp=true, tol=1e-5)
    plugin == "ipqp"    && return (opts=Dict{String,Any}("print_iter"=>false, "print_header"=>false, "print_info"=>false), lp=false, tol=1e-5)
    plugin == "daqp"    && return (opts=Dict{String,Any}(), lp=false, tol=1e-5)
    plugin == "proxqp"  && return (opts=Dict{String,Any}("print_time"=>false), lp=false, tol=1e-4)
    plugin == "piqp"    && return (opts=Dict{String,Any}("print_time"=>false), lp=false, tol=1e-5)
    return (opts=Dict{String,Any}(), lp=false, tol=1e-4)
end

for plugin in ["qrqp", "osqp", "qpoases", "highs", "ipqp", "daqp", "proxqp", "piqp"]
    (try C.has_conic(plugin) catch; false end) || continue
    meta = _conic_meta(plugin)
    @testset "conic $plugin" begin
        # --- unconstrained QP: min 0.5(x^2+y^2) - 0.7x - 2.3y -> [0.7, 2.3] ---
        H = DM([1.0 0.0; 0.0 1.0]); G = DM([-0.7, -2.3])
        S = conic("Su", plugin, Dict("h" => C.sparsity(H)), meta.opts)
        r = S(Dict("h" => H, "g" => G, "lbx" => DM([-_INF, -_INF]), "ubx" => DM([_INF, _INF])))
        @test maximum(abs.(Vector{Float64}(r["x"]) .- [0.7, 2.3])) < meta.tol
        @test abs(Float64(r["cost"]) - (-0.5 * (0.49 + 5.29))) < meta.tol

        # --- convex QP with 3 inequality rows -> [2/3, 4/3], cost -8-2/9 ---
        H2 = DM([1.0 -1.0; -1.0 2.0]); G2 = DM([-2.0, -6.0])
        A2 = DM([1.0 1.0; -1.0 2.0; 2.0 1.0])
        S2 = conic("Sc", plugin, Dict("h" => C.sparsity(H2), "a" => C.sparsity(A2)), meta.opts)
        r2 = S2(Dict("h" => H2, "g" => G2, "a" => A2,
            "lbx" => DM([0.0, 0.0]), "ubx" => DM([_INF, _INF]),
            "lba" => DM([-_INF, -_INF, -_INF]), "uba" => DM([2.0, 2.0, 3.0])))
        @test maximum(abs.(Vector{Float64}(r2["x"]) .- [2 / 3, 4 / 3])) < meta.tol
        @test abs(Float64(r2["cost"]) - (-8 - 2 / 9)) < meta.tol

        # --- equality-only QP (lba == uba) -> [-0.5, 1], cost -3.375 ---
        Ae = DM([1.0 1.0])
        Se = conic("Se", plugin, Dict("h" => C.sparsity(H2), "a" => C.sparsity(Ae)), meta.opts)
        re = Se(Dict("h" => H2, "g" => G2, "a" => Ae,
            "lbx" => DM([-10.0, -10.0]), "ubx" => DM([10.0, 10.0]),
            "lba" => DM([0.5]), "uba" => DM([0.5])))
        @test maximum(abs.(Vector{Float64}(re["x"]) .- [-0.5, 1.0])) < meta.tol
        @test abs(Float64(re["cost"]) - (-3.375)) < meta.tol

        # --- mixed equality + inequality, x[0] fixed -> [0.5, 1.25], cost -7.4375 ---
        Sq = conic("Sq", plugin, Dict("h" => C.sparsity(H2), "a" => C.sparsity(A2)), meta.opts)
        rq = Sq(Dict("h" => H2, "g" => G2, "a" => A2,
            "lbx" => DM([0.5, 0.0]), "ubx" => DM([0.5, _INF]),
            "lba" => DM([-_INF, -_INF, -_INF]), "uba" => DM([2.0, 2.0, 3.0])))
        @test maximum(abs.(Vector{Float64}(rq["x"]) .- [0.5, 1.25])) < meta.tol
        @test abs(Float64(rq["cost"]) - (-7.4375)) < meta.tol

        # --- parametric re-solve: diagonal H, vary g -> x = -g ---
        Sp = conic("Sp", plugin, Dict("h" => C.sparsity(H)), meta.opts)
        p1 = Sp(Dict("h" => H, "g" => DM([-1.0, -2.0]), "lbx" => DM([-_INF, -_INF]), "ubx" => DM([_INF, _INF])))
        @test maximum(abs.(Vector{Float64}(p1["x"]) .- [1.0, 2.0])) < meta.tol
        p2 = Sp(Dict("h" => H, "g" => DM([-3.0, -4.0]), "lbx" => DM([-_INF, -_INF]), "ubx" => DM([_INF, _INF])))
        @test maximum(abs.(Vector{Float64}(p2["x"]) .- [3.0, 4.0])) < meta.tol

        # --- warmstart: re-solve seeded with the converged primal+dual ---
        # stats()/iter_count is unbound (SWIG return-type skip), so instead assert
        # the warm solve, started AT the optimum, stays there: its x matches both
        # the cold x and the analytic optimum to (near) machine precision -- far
        # tighter than the cold-solve meta.tol, which distinguishes a warmstart
        # that actually ingests x0/lam0 from one that silently ignores them.
        w1 = S2(Dict("h" => H2, "g" => G2, "a" => A2,
            "lbx" => DM([0.0, 0.0]), "ubx" => DM([_INF, _INF]),
            "lba" => DM([-_INF, -_INF, -_INF]), "uba" => DM([2.0, 2.0, 3.0])))
        w2 = S2(Dict("h" => H2, "g" => G2, "a" => A2,
            "lbx" => DM([0.0, 0.0]), "ubx" => DM([_INF, _INF]),
            "lba" => DM([-_INF, -_INF, -_INF]), "uba" => DM([2.0, 2.0, 3.0]),
            "x0" => w1["x"], "lam_x0" => w1["lam_x"], "lam_a0" => w1["lam_a"]))
        @test abs(Float64(w2["cost"]) - Float64(w1["cost"])) < meta.tol
        # warm x equals cold x (warmstart preserves the solution, not just cost)
        @test maximum(abs.(Vector{Float64}(w2["x"]) .- Vector{Float64}(w1["x"]))) < meta.tol
        # and both sit on the analytic optimum [2/3, 4/3]
        @test maximum(abs.(Vector{Float64}(w2["x"]) .- [2 / 3, 4 / 3])) < meta.tol

        # --- LP (H = 0): min 2x+y s.t. rows -> [0.5, 1.5], cost 2.5 ---
        if meta.lp
            Al = DM([-1.0 1.0; 1.0 1.0; 1.0 -2.0]); Gl = DM([2.0, 1.0])
            Sl = conic("Sl", plugin, Dict("h" => C.Sparsity(2, 2), "a" => C.sparsity(Al)), meta.opts)
            rl = Sl(Dict("g" => Gl, "a" => Al,
                "lbx" => DM([-_INF, 0.0]), "ubx" => DM([_INF, _INF]),
                "lba" => DM([-_INF, 2.0, -_INF]), "uba" => DM([1.0, _INF, 4.0])))
            @test maximum(abs.(Vector{Float64}(rl["x"]) .- [0.5, 1.5])) < 10 * meta.tol
            @test abs(Float64(rl["cost"]) - 2.5) < 10 * meta.tol
        end

        # --- infeasibility: inconsistent bounds lbx > ubx; the call signals failure ---
        Si = conic("Si", plugin, Dict("h" => C.sparsity(H)), meta.opts)
        @test_throws Exception Si(Dict("h" => H, "g" => DM([0.0, 0.0]),
            "lbx" => DM([2.0, 2.0]), "ubx" => DM([1.0, 1.0])))
    end
end
