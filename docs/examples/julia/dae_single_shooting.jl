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
# Compact direct single shooting for a DAE optimal control problem:
#   minimize     integral_0^10 x0^2 + x1^2 + u^2 dt
#   subject to   dot(x0) == z*x0 - x1 + u
#                dot(x1) == x0
#                      0 == x1^2 + z - 1
#                x0(0)==0, x1(0)==1, x0(10)==0, x1(10)==0, -0.75 <= u <= 1
#
# "collocation" stands in for the python "idas" integrator (no sundials here).
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

# Declare variables
x = ca.SX.sym("x", 2)   # Differential states
z = ca.SX.sym("z")      # Algebraic variable
u = ca.SX.sym("u")      # Control

# Differential equation
f_x = ca.vertcat(z * x[1] - x[2] + u, x[1])

# Algebraic equation
f_z = x[2]^2 + z - 1

# Lagrange cost term (quadrature)
f_q = x[1]^2 + x[2]^2 + u^2

# Create an integrator (interval length 0.5s)
dae = Dict("x" => x, "z" => z, "p" => u, "ode" => f_x, "alg" => f_z, "quad" => f_q)
I = ca.integrator("I", "collocation", dae, 0.0, 0.5)

# All controls
U = ca.MX.sym("U", 20)

function main()
    # Construct graph of integrator calls
    X = ca.MX(ca.DM([0.0, 1.0]))
    J = ca.MX(ca.DM(0.0))
    for k in 1:20
        Ik = I(x0 = X, p = U[k])
        X = Ik["xf"]
        J = J + Ik["qf"]      # Sum up quadratures
    end

    # Allocate an NLP solver
    nlp = Dict("x" => U, "f" => J, "g" => X)
    solver = ca.nlpsol("solver", "ipopt", nlp)

    # Pass bounds, initial guess and solve NLP
    sol = solver(lbx = -0.75,   # Lower variable bound
                 ubx = 1.0,     # Upper variable bound
                 lbg = 0.0,     # Lower constraint bound
                 ubg = 0.0,     # Upper constraint bound
                 x0 = 0.0)      # Initial guess

    println("-----")
    println("objective at solution = ", sol["f"])
    println("u_opt = ", sol["x"])
    println("terminal state x(10) = ", sol["g"])
end

main()
