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
# Direct collocation transcription of a Van der Pol OCP, building the
# collocation coefficients by hand (no integrator plugin) -> fully portable.
# numpy poly1d/polyder/polyint are replaced by small array-based helpers
# (coefficients highest-degree-first, like numpy.poly1d).
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

# Polynomial helpers (coeffs highest-degree-first, like numpy.poly1d)
function polymul(a, b)
    r = zeros(length(a) + length(b) - 1)
    for i in 1:length(a), j in 1:length(b)
        r[i + j - 1] += a[i] * b[j]
    end
    r
end
polyval(p, x) = foldl((r, c) -> r * x + c, p; init = 0.0)
polyder(p) = (n = length(p) - 1; n == 0 ? [0.0] : [p[i] * (n - i + 1) for i in 1:n])
polyint(p) = vcat([p[i] / (length(p) - i + 1) for i in 1:length(p)], 0.0)

function main()
    # Degree of interpolating polynomial
    d = 3

    # Collocation points (prepend 0)
    tau_root = [0.0; ca.collocation_points(d, "legendre")]

    # Coefficients of the collocation, continuity and quadrature equations
    Cc = zeros(d + 1, d + 1)
    D = zeros(d + 1)
    B = zeros(d + 1)

    # Construct polynomial basis
    for j in 1:(d + 1)
        p = [1.0]
        for r in 1:(d + 1)
            if r != j
                p = polymul(p, [1.0, -tau_root[r]]) ./ (tau_root[j] - tau_root[r])
            end
        end
        # Continuity coefficient
        D[j] = polyval(p, 1.0)
        # Collocation coefficients
        pder = polyder(p)
        for r in 1:(d + 1)
            Cc[j, r] = polyval(pder, tau_root[r])
        end
        # Quadrature coefficient
        B[j] = polyval(polyint(p), 1.0)
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
    N = 20      # number of control intervals
    h = T / N

    # Start with an empty NLP
    w = ca.MX[]; w0 = Float64[]; lbw = Float64[]; ubw = Float64[]
    J = ca.MX(ca.DM(0.0))
    g = ca.MX[]; lbg = Float64[]; ubg = Float64[]

    # For plotting/logging x and u given w
    x_plot = ca.MX[]; u_plot = ca.MX[]

    # "Lift" initial conditions
    Xk = ca.MX.sym("X0", 2)
    push!(w, Xk); append!(lbw, [0, 1]); append!(ubw, [0, 1]); append!(w0, [0, 1])
    push!(x_plot, Xk)

    # Formulate the NLP
    for k in 0:(N - 1)
        # New NLP variable for the control
        Uk = ca.MX.sym("U_$k")
        push!(w, Uk); push!(lbw, -1); push!(ubw, 1); push!(w0, 0)
        push!(u_plot, Uk)

        # State at collocation points
        Xc = ca.MX[]
        for j in 1:d
            Xkj = ca.MX.sym("X_$(k)_$(j)", 2)
            push!(Xc, Xkj); push!(w, Xkj)
            append!(lbw, [-0.25, -Inf]); append!(ubw, [Inf, Inf]); append!(w0, [0, 0])
        end

        # Loop over collocation points
        Xk_end = D[1] * Xk
        for j in 1:d
            # Expression for the state derivative at the collocation point
            xp = Cc[1, j + 1] * Xk
            for r in 1:d
                xp = xp + Cc[r + 1, j + 1] * Xc[r]
            end

            # Append collocation equations
            fj, qj = f(Xc[j], Uk)
            push!(g, h * fj - xp); append!(lbg, [0, 0]); append!(ubg, [0, 0])

            # Add contribution to the end state
            Xk_end = Xk_end + D[j + 1] * Xc[j]

            # Add contribution to quadrature function
            J = J + B[j + 1] * qj * h
        end

        # New NLP variable for state at end of interval
        Xk = ca.MX.sym("X_$(k + 1)", 2)
        push!(w, Xk)
        append!(lbw, [-0.25, -Inf]); append!(ubw, [Inf, Inf]); append!(w0, [0, 0])
        push!(x_plot, Xk)

        # Add equality constraint
        push!(g, Xk_end - Xk); append!(lbg, [0, 0]); append!(ubg, [0, 0])
    end

    # Create an NLP solver
    prob = Dict("f" => J, "x" => ca.vertcat(w), "g" => ca.vertcat(g))
    solver = ca.nlpsol("solver", "ipopt", prob)

    # Function to get x and u trajectories from w
    W = ca.vertcat(w)
    trajectories = ca.Function("trajectories", [W], [ca.horzcat(x_plot), ca.horzcat(u_plot)],
                              ["w"], ["x", "u"])

    # Solve the NLP
    sol = solver(x0 = w0, lbx = lbw, ubx = ubw, lbg = lbg, ubg = ubg)
    x_opt, u_opt = trajectories(sol["x"])

    println("-----")
    println("objective at solution = ", sol["f"])
    println("x1 trajectory = ", x_opt[1, :])
    println("x2 trajectory = ", x_opt[2, :])
    println("u trajectory = ", u_opt)
end

main()
