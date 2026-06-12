# linsol.jl -- capability port of test/javascript/linsol.js (+ deeper).
# For every loaded linear-solver plugin: factorize and solve A x = b,
# transpose-solve, multiple RHS, symbolic sparsity reuse, MX-level solve,
# residual checks, and a singular-matrix error path.

const _LIN_OD = Dict{String,Any}()

# Per-solver requirements (symmetry/posdef/precision) mirroring linsol.js.
const LIN_SOLVERS = [
    ("qr",         (sym=false, iter=false, tol=1e-8)),
    ("ldl",        (sym=true,  iter=false, tol=1e-8)),
    ("lsqr",       (sym=false, iter=true,  tol=1e-4)),
    ("symbolicqr", (sym=false, iter=false, tol=1e-8)),
    ("csparse",    (sym=false, iter=false, tol=1e-8)),
    ("lapacklu",   (sym=false, iter=false, tol=1e-8)),
    ("lapackqr",   (sym=false, iter=false, tol=1e-8)),
    ("ldl",        (sym=true,  iter=false, tol=1e-8)),
]

# Reference 2x2 solve via Cramer's rule on a row-major [[a,b],[c,d]].
_solve2(A, b) = (det = A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1];
    [(A[2, 2] * b[1] - A[1, 2] * b[2]) / det,
     (-A[2, 1] * b[1] + A[1, 1] * b[2]) / det])

_makeA2(sym) = sym ? DM([4.0 1.0; 1.0 3.0]) : DM([3.0 7.0; 1.0 2.0])
_Arow2(sym) = sym ? [4.0 1.0; 1.0 3.0] : [3.0 7.0; 1.0 2.0]

for (plugin, req) in LIN_SOLVERS
    (try C.has_linsol(plugin) catch; false end) || continue
    @testset "linsol $plugin" begin
        # --- basic dense 2x2 via the Linsol class ---
        A = _makeA2(req.sym); b = DM([1.0, 0.5])
        ls = C.Linsol("S", plugin, C.sparsity(A), _LIN_OD)
        x = C.solve(ls, A, b, false)
        ref = _solve2(_Arow2(req.sym), [1.0, 0.5])
        @test maximum(abs.(Vector{Float64}(x) .- ref)) < req.tol

        # --- transpose solve A' x = b ---
        xt = C.solve(ls, A, b, true)
        reft = _solve2(permutedims(_Arow2(req.sym)), [1.0, 0.5])
        @test maximum(abs.(Vector{Float64}(xt) .- reft)) < req.tol

        # --- explicit sfact/nfact split, then solve ---
        ls2 = C.Linsol("S", plugin, C.sparsity(A), _LIN_OD)
        C.sfact(ls2, A); C.nfact(ls2, A)
        x2 = C.solve(ls2, A, b, false)
        @test maximum(abs.(Vector{Float64}(x2) .- ref)) < req.tol

        # --- multiple RHS: A X = B, check residual ---
        B = DM([1.0 0.3; 0.5 0.7])
        X = C.solve(ls, A, B, false)
        @test Base.size(X) == (2, 2)
        resid = Matrix{Float64}(C.minus(C.mtimes(A, X), B))
        @test maximum(abs.(resid)) < req.tol

        # --- symbolic sparsity reuse: refactor a different A, same pattern ---
        A3 = DM(req.sym ? [8.0 2.0; 2.0 6.0] : [6.0 14.0; 2.0 4.0])
        b3 = DM([2.0, 1.0])
        x3 = C.solve(ls, A3, b3, false)
        ref3 = _solve2(req.sym ? [8.0 2.0; 2.0 6.0] : [6.0 14.0; 2.0 4.0], [2.0, 1.0])
        @test maximum(abs.(Vector{Float64}(x3) .- ref3)) < req.tol

        # --- free-function solve(A, b, plugin, opts) DM path ---
        xf = C.solve(A, b, plugin, _LIN_OD)
        @test maximum(abs.(Vector{Float64}(xf) .- ref)) < req.tol

        # --- MX-level solve wrapped in a Function, evaluated numerically ---
        Am = sym(MX, "A", 2, 2); bm = sym(MX, "b", 2)
        xm = C.solve(Am, bm, plugin, _LIN_OD)
        fm = C.Function("f", [Am, bm], [xm])
        outm = fm(A, b)
        @test maximum(abs.(Vector{Float64}(outm) .- ref)) < req.tol

        # --- identity & diagonal sanity ---
        I3 = DM([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
        lsI = C.Linsol("I", plugin, C.sparsity(I3), _LIN_OD)
        xi = C.solve(lsI, I3, DM([3.0, 5.0, 7.0]), false)
        @test maximum(abs.(Vector{Float64}(xi) .- [3.0, 5.0, 7.0])) < req.tol
        Dg = DM([2.0 0.0 0.0; 0.0 3.0 0.0; 0.0 0.0 4.0])
        lsD = C.Linsol("D", plugin, C.sparsity(Dg), _LIN_OD)
        xd = C.solve(lsD, Dg, DM([4.0, 9.0, 16.0]), false)
        @test maximum(abs.(Vector{Float64}(xd) .- [2.0, 3.0, 4.0])) < req.tol

        # --- 4x4 well-conditioned residual ---
        A4 = req.sym ? DM([5.0 1.0 0.0 0.0; 1.0 4.0 1.0 0.0; 0.0 1.0 3.0 1.0; 0.0 0.0 1.0 2.0]) :
                       DM([3.0 1.0 0.0 0.0; 1.0 4.0 1.0 0.0; 0.0 1.0 5.0 1.0; 0.0 0.0 1.0 6.0])
        b4 = DM([1.0, 2.0, 3.0, 4.0])
        ls4 = C.Linsol("S4", plugin, C.sparsity(A4), _LIN_OD)
        x4 = C.solve(ls4, A4, b4, false)
        r4 = Vector{Float64}(C.minus(C.mtimes(A4, x4), b4))
        @test maximum(abs.(r4)) < req.tol

        # --- singular matrix: direct LU/QR factorizers detect rank deficiency ---
        if plugin in ("qr", "csparse", "lapacklu")
            Asing = DM([1.0 2.0; 2.0 4.0])
            lss = C.Linsol("Sing", plugin, C.sparsity(Asing), _LIN_OD)
            C.sfact(lss, Asing)
            @test_throws Exception C.nfact(lss, Asing)
        end
    end
end
