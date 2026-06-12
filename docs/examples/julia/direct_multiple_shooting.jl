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
# Direct multiple shooting of the Van der Pol OCP, with a fixed-step RK4
# integrator built by hand and lifted interval states.
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

T = 10.0    # Time horizon
N = 20      # number of control intervals

# Declare model variables
x1 = ca.MX.sym("x1")
x2 = ca.MX.sym("x2")
x = ca.vertcat(x1, x2)
u = ca.MX.sym("u")

# Model equations
xdot = ca.vertcat((1 - x2^2) * x1 - x2 + u, x1)

# Objective term
L = x1^2 + x2^2 + u^2

# Formulate discrete time dynamics: fixed step Runge-Kutta 4 integrator
M = 4              # RK4 steps per interval
DT = T / N / M
f = ca.Function("f", [x, u], [xdot, L])

function rk4_step(F)
    X0 = ca.MX.sym("X0", 2)
    U = ca.MX.sym("U")
    X = X0
    Q = ca.MX(ca.DM(0.0))
    for j in 1:M
        k1, k1q = F(X, U)
        k2, k2q = F(X + DT / 2 * k1, U)
        k3, k3q = F(X + DT / 2 * k2, U)
        k4, k4q = F(X + DT * k3, U)
        X = X + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        Q = Q + DT / 6 * (k1q + 2 * k2q + 2 * k3q + k4q)
    end
    ca.Function("F", [X0, U], [X, Q], ["x0", "p"], ["xf", "qf"])
end

function main()
    F = rk4_step(f)

    # Evaluate at a test point
    Fk = F(x0 = [0.2, 0.3], p = 0.4)
    println("test xf = ", Fk["xf"])
    println("test qf = ", Fk["qf"])

    # Start with an empty NLP
    w = ca.MX[]; w0 = Float64[]; lbw = Float64[]; ubw = Float64[]
    J = ca.MX(ca.DM(0.0))
    g = ca.MX[]; lbg = Float64[]; ubg = Float64[]

    # "Lift" initial conditions
    Xk = ca.MX.sym("X0", 2)
    push!(w, Xk); append!(lbw, [0, 1]); append!(ubw, [0, 1]); append!(w0, [0, 1])

    # Formulate the NLP
    for k in 0:N - 1
        # New NLP variable for the control
        Uk = ca.MX.sym("U_$k")
        push!(w, Uk); push!(lbw, -1); push!(ubw, 1); push!(w0, 0)

        # Integrate till the end of the interval
        Fk = F(x0 = Xk, p = Uk)
        Xk_end = Fk["xf"]
        J = J + Fk["qf"]

        # New NLP variable for state at end of interval
        Xk = ca.MX.sym("X_$(k + 1)", 2)
        push!(w, Xk)
        append!(lbw, [-0.25, -Inf]); append!(ubw, [Inf, Inf]); append!(w0, [0, 0])

        # Add equality constraint
        push!(g, Xk_end - Xk); append!(lbg, [0, 0]); append!(ubg, [0, 0])
    end

    # Create an NLP solver
    prob = Dict("f" => J, "x" => ca.vertcat(w), "g" => ca.vertcat(g))
    solver = ca.nlpsol("solver", "ipopt", prob)

    # Solve the NLP
    sol = solver(x0 = w0, lbx = lbw, ubx = ubw, lbg = lbg, ubg = ubg)
    w_opt = sol["x"]

    println("-----")
    println("objective at solution = ", sol["f"])
    println("x1_opt = ", w_opt[1:3:end])
    println("x2_opt = ", w_opt[2:3:end])
    println("u_opt  = ", w_opt[3:3:end])
end

main()
