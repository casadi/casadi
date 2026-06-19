# rootfinder.jl -- capability port of test/javascript/rootfinder.js (+ deeper).
# For every loaded rootfinder plugin: scalar and 2-D nonlinear systems,
# parametric residuals, multimodal convergence, and an embedded (MX-level)
# rootfinder.  Newton uses the analytic Jacobian of the residual.

const _RF_OD = Dict{String,Any}()

# Per-plugin construction options (some need a sub-solver / quiet flags).
function _rf_opts(plugin)
    plugin == "nlpsol"      && return Dict{String,Any}("nlpsol" => "ipopt",
        "nlpsol_options" => Dict{String,Any}("print_time" => false,
            "ipopt" => Dict{String,Any}("print_level" => 0, "sb" => "yes")))
    plugin == "newton"      && return Dict{String,Any}("print_iteration" => false)
    plugin == "kinsol"      && return Dict{String,Any}("print_level" => 0)
    return _RF_OD
end

for plugin in ["newton", "fast_newton", "kinsol", "nlpsol"]
    (try C.has_rootfinder(plugin) catch; false end) || continue
    plugin == "nlpsol" && !HAS_IPOPT && continue   # rootfinder:nlpsol uses ipopt inside
    @testset "rootfinder $plugin" begin
        # --- sin(x) = 0, x0 = 6 -> 2*pi ---
        x = sym(SX, "x")
        f = C.Function("f", [x], [sin(x)])
        F = C.rootfinder("F", plugin, f, _rf_opts(plugin))
        @test abs(Float64(F(DM(6.0))) - 2pi) < 1e-6

        # --- multimodal: same solver, two start points ---
        @test abs(Float64(F(DM(1.0))) - 0.0) < 1e-5
        @test abs(Float64(F(DM(4.0))) - pi) < 1e-5

        # --- x^2 - 2 = 0, x0 = 1.5 -> sqrt(2) ---
        g = C.Function("g", [x], [C.times(x, x) - 2.0])
        G = C.rootfinder("G", plugin, g, _rf_opts(plugin))
        @test abs(Float64(G(DM(1.5))) - sqrt(2)) < 1e-6

        # --- cubic x^3 - x - 1 = 0 -> 1.324717957... ---
        x3 = C.times(C.times(x, x), x)
        h = C.Function("h", [x], [x3 - x - 1.0])
        H = C.rootfinder("H", plugin, h, _rf_opts(plugin))
        @test abs(Float64(H(DM(1.5))) - 1.3247179572447460) < 1e-5

        # --- parametric: residual x_par - asin(y) = 0 -> y = sin(x_par) ---
        y = sym(SX, "y"); xp = sym(SX, "xp")
        fp = C.Function("fp", [y, xp], [xp - C.asin(y)])
        Fp = C.rootfinder("Fp", plugin, fp, _rf_opts(plugin))
        @test abs(Float64(Fp(DM(0.0), DM(0.2))) - sin(0.2)) < 1e-6

        # --- 2-D system: x^2 + y^2 - 2 = 0, x - y = 0 -> (1, 1) ---
        xv = sym(SX, "x"); yv = sym(SX, "y")
        xy = vertcat(xv, yv)
        res = vertcat(C.times(xv, xv) + C.times(yv, yv) - 2.0, xv - yv)
        f2 = C.Function("f2", [xy], [res])
        F2 = C.rootfinder("F2", plugin, f2, _rf_opts(plugin))
        r2 = Vector{Float64}(F2(DM([0.5, 0.5])))
        @test maximum(abs.(r2 .- [1.0, 1.0])) < 1e-5

        # --- 2-D linear system with A, b as parameters -> A \ b ---
        xpar = sym(SX, "x", 2); A = sym(SX, "A", 2, 2); b = sym(SX, "b", 2)
        flin = C.Function("flin", [xpar, A, b], [C.mtimes(A, xpar) - b])
        Flin = C.rootfinder("Flin", plugin, flin, _rf_opts(plugin))
        Av = DM([1.0 2.0; 3.0 2.1]); bv = DM([0.7, 0.6])
        rlin = Vector{Float64}(Flin(DM([0.0, 0.0]), Av, bv))
        det = 1 * 2.1 - 2 * 3
        @test maximum(abs.(rlin .- [(2.1 * 0.7 - 2 * 0.6) / det, (0.6 - 3 * 0.7) / det])) < 1e-6

        # --- embedded in another Function, evaluated at MX level ---
        X = sym(MX, "X")
        R = Fp(C.MX(0.0), X)
        trial = C.Function("trial", [X], [R])
        @test abs(Float64(trial(DM(0.2))) - sin(0.2)) < 1e-6
    end
end
