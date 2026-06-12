#
#     MIT No Attribution
#
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# Direct-collocation Van der Pol OCP solved with an SQP method (qrqp QP
# sub-solver) for accurate multipliers, then NLP-solution sensitivities w.r.t. a
# perturbation parameter P computed three ways:
#   - high-level: factory for the Hessian of f w.r.t. p,
#   - low-level forward AD (solver.forward) vs finite differences,
#   - low-level reverse AD (solver.reverse).
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

# Degree of interpolating polynomial
d = 3
# Get collocation points
tau_root = vcat(0.0, ca.collocation_points(d, "legendre"))

# Collocation (Coll), continuity (Cont) and quadrature (Quad) coefficients,
# obtained by evaluating the Lagrange basis (and its derivative / integral).
Coll = zeros(d + 1, d + 1)
Cont = zeros(d + 1)
Quad = zeros(d + 1)
tau = ca.SX.sym("tau")
for j in 1:d + 1
    L = ca.SX(1.0)
    for r in 1:d + 1
        if r != j
            L = L * (tau - tau_root[r]) / (tau_root[j] - tau_root[r])
        end
    end
    Cont[j] = Float64(ca.Function("l", [tau], [L])(1.0))
    dL = ca.tangent(L, tau)
    tfcn = ca.Function("t", [tau], [dL])
    for r in 1:d + 1
        Coll[j, r] = Float64(tfcn(tau_root[r]))
    end
    # integral of the basis over [0,1] via a fine Simpson-free quadrature:
    # build an antiderivative symbolically is awkward, so integrate the
    # polynomial numerically with a dense rk integrator.
    s = ca.SX.sym("s")
    iquad = ca.integrator("iq", "rk", Dict("x" => s, "t" => tau, "ode" => L),
                       0.0, 1.0, Dict("number_of_finite_elements" => 50))
    Quad[j] = Float64(iquad(x0 = 0.0)["xf"])
end

# Time horizon
T = 10.0

# Declare model variables
x1 = ca.SX.sym("x1")
x2 = ca.SX.sym("x2")
x = ca.vertcat(x1, x2)
u = ca.SX.sym("u")

# Model equations
xdot = ca.vertcat((1 - x2^2) * x1 - x2 + u, x1)

# Objective term
L = x1^2 + x2^2 + u^2

# Continuous time dynamics
f = ca.Function("f", [x, u], [xdot, L], ["x", "u"], ["xdot", "L"])

# Control discretization
N = 20    # number of control intervals
h = T / N

function build(P)
    # Start with an empty NLP
    w = ca.MX[]; w0 = Float64[]; lbw = Float64[]; ubw = Float64[]
    J = ca.MX(ca.DM(0.0))
    g = ca.MX[]; lbg = Float64[]; ubg = Float64[]
    x_plot = ca.MX[]; u_plot = ca.MX[]

    # "Lift" initial conditions
    Xk = ca.MX.sym("X0", 2)
    push!(w, Xk); append!(lbw, [0, 1]); append!(ubw, [0, 1]); append!(w0, [0, 1])
    push!(x_plot, Xk)

    # Perturb with P
    Xk = Xk + P

    for k in 0:N - 1
        # New NLP variable for the control
        Uk = ca.MX.sym("U_$k")
        push!(w, Uk); push!(lbw, -1); push!(ubw, 0.85); push!(w0, 0)
        push!(u_plot, Uk)

        # State at collocation points
        Xc = ca.MX[]
        for j in 1:d
            Xkj = ca.MX.sym("X_$(k)_$(j)", 2)
            push!(Xc, Xkj); push!(w, Xkj)
            append!(lbw, [-0.25, -Inf]); append!(ubw, [Inf, Inf]); append!(w0, [0, 0])
        end

        # Loop over collocation points
        Xk_end = Cont[1] * Xk
        for j in 2:d + 1
            xp = Coll[1, j] * Xk
            for r in 1:d
                xp = xp + Coll[r + 1, j] * Xc[r]
            end
            fj, qj = f(Xc[j - 1], Uk)
            push!(g, h * fj - xp); append!(lbg, [0, 0]); append!(ubg, [0, 0])
            Xk_end = Xk_end + Cont[j] * Xc[j - 1]
            J = J + Quad[j] * qj * h
        end

        # New NLP variable for the state at the end of the interval
        Xk = ca.MX.sym("X_$(k + 1)", 2)
        push!(w, Xk)
        append!(lbw, [-0.25, -Inf]); append!(ubw, [Inf, Inf]); append!(w0, [0, 0])
        push!(x_plot, Xk)

        # Continuity (equality) constraint
        push!(g, Xk_end - Xk); append!(lbg, [0, 0]); append!(ubg, [0, 0])
    end

    (ca.vertcat(w), J, ca.vertcat(g), w0, lbw, ubw, lbg, ubg)
end

function main()
    P = ca.MX.sym("P", 2)
    w, J, g, w0, lbw, ubw, lbg, ubg = build(P)

    # NLP with the perturbation parameter
    prob = Dict("f" => J, "x" => w, "g" => g, "p" => P)

    # SQP with active-set QP for accurate multipliers
    solver = ca.nlpsol("solver", "sqpmethod", prob, Dict(
        "qpsol" => "qrqp",
        "qpsol_options" => Dict{String,Any}("print_iter" => false, "error_on_fail" => false),
        "print_time" => false))

    # Solve the NLP at the nominal parameter value
    sol = solver(x0 = w0, lbx = lbw, ubx = ubw, lbg = lbg, ubg = ubg, p = 0)

    # High-level approach: Hessian of optimal f w.r.t. p via factory
    hsolver = ca.factory(solver, "h", ca.name_in(solver), ["hess:f:p:p"])
    println("hsolver generated")
    hsol = hsolver(x0 = sol["x"], lam_x0 = sol["lam_x"], lam_g0 = sol["lam_g"],
                   lbx = lbw, ubx = ubw, lbg = lbg, ubg = ubg, p = 0)
    println("Hessian of f w.r.t. p:")
    println(hsol["hess_f_p_p"])

    # Low-level: forward-mode directional derivatives, two directions at once
    nfwd = 2
    fwd_solver = solver.forward(nfwd)
    println("fwd_solver generated")

    # Seeds: perturb P in two directions (its first and second nonzero)
    spx = ca.sparsity(sol["x"]); spg = ca.sparsity(sol["g"]); spP = ca.sparsity(P)
    fwd_p = [ca.zeros(ca.GenDM, spP) for _ in 1:nfwd]
    fwd_p[1][1] = 1
    fwd_p[2][2] = 1
    # Typed horzcat keeps the all-DM seed stacks numeric.
    z(sp) = ca.horzcat([ca.zeros(ca.GenDM, sp) for _ in 1:nfwd])

    sol_fwd = fwd_solver(out_x = sol["x"], out_lam_g = sol["lam_g"], out_lam_x = sol["lam_x"],
                         out_f = sol["f"], out_g = sol["g"], lbx = lbw, ubx = ubw, lbg = lbg, ubg = ubg,
                         fwd_lbx = z(spx), fwd_ubx = z(spx), fwd_lbg = z(spg), fwd_ubg = z(spg),
                         p = 0, fwd_p = ca.horzcat(fwd_p))

    # Same thing via finite differences
    hfd = 1e-3
    pert = [solver(x0 = sol["x"], lam_g0 = sol["lam_g"], lam_x0 = sol["lam_x"],
                   lbx = lbw, ubx = ubw, lbg = lbg, ubg = ubg, p = 0 + hfd * fwd_p[dd])
            for dd in 1:nfwd]

    println("==========")
    println("Checking f")
    println("finite differences")
    for dd in 1:nfwd
        println((pert[dd]["f"] - sol["f"]) / hfd)
    end
    println("AD fwd")
    M = sol_fwd["fwd_f"]
    for dd in 1:nfwd
        println(M[1, dd])
    end

    # Second-order finite differences (perturb the opposite way too)
    pert2 = [solver(x0 = sol["x"], lam_g0 = sol["lam_g"], lam_x0 = sol["lam_x"],
                    lbx = lbw, ubx = ubw, lbg = lbg, ubg = ubg, p = 0 - hfd * fwd_p[dd])
             for dd in 1:nfwd]
    println("finite differences, second order: f")
    for dd in 1:nfwd
        println((pert[dd]["f"] - 2 * sol["f"] + pert2[dd]["f"]) / (hfd * hfd))
    end

    # Low-level: reverse-mode AD for the NLP solver object
    nadj = 1
    adj_solver = solver.reverse(nadj)   # `ca.reverse` would collide with Base.reverse
    println("adj_solver generated")
    spf = ca.sparsity(sol["f"])
    adj_f = ca.zeros(ca.GenDM, spf); adj_f[1] = 1
    za(sp) = ca.zeros(ca.GenDM, sp)

    sol_adj = adj_solver(out_x = sol["x"], out_lam_g = sol["lam_g"], out_lam_x = sol["lam_x"],
                         out_f = sol["f"], out_g = sol["g"], lbx = lbw, ubx = ubw, lbg = lbg, ubg = ubg,
                         adj_f = adj_f, adj_g = za(spg), p = 0, adj_x = za(spx))

    println("==========")
    println("Checking p")
    println("Reverse mode AD")
    println(sol_adj["adj_p"])
end

main()
