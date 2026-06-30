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
# Direct multiple shooting for the same DAE optimal control problem as
# dae_single_shooting. "collocation" stands in for the python "idas" integrator.
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

# Declare variables
x0 = ca.SX.sym("x0")
x1 = ca.SX.sym("x1")
x = ca.vertcat(x0, x1)   # Differential states
z = ca.SX.sym("z")      # Algebraic variable
u = ca.SX.sym("u")      # Control

# Differential equation
f_x = ca.vertcat(z * x0 - x1 + u, x0)

# Algebraic equation
f_z = x1^2 + z - 1

# Lagrange cost term (quadrature)
f_q = x0^2 + x1^2 + u^2

# Create an integrator (interval length 0.5s)
dae = Dict("x" => x, "z" => z, "p" => u, "ode" => f_x, "alg" => f_z, "quad" => f_q)
I = ca.integrator("I", "collocation", dae, 0.0, 0.5)

# Number of intervals
nk = 20

function main()
    # Start with an empty NLP
    w = ca.MX[]              # List of variables
    lbw = Float64[]       # Lower bounds on w
    ubw = Float64[]       # Upper bounds on w
    G = ca.MX[]              # Constraints
    J = ca.MX(ca.DM(0.0))   # Cost function

    # Initial conditions
    Xk = ca.MX.sym("X0", 2)
    push!(w, Xk)
    append!(lbw, [0, 1]); append!(ubw, [0, 1])

    # Loop over all intervals
    for k in 0:nk - 1
        # Local control
        Uk = ca.MX.sym("U$k")
        push!(w, Uk)
        push!(lbw, -0.75); push!(ubw, 1.00)

        # Call integrator function
        Ik = I(x0 = Xk, p = Uk)
        Xk = Ik["xf"]
        J = J + Ik["qf"]   # Sum quadratures

        # "Lift" the variable
        X_prev = Xk
        Xk = ca.MX.sym("X$(k + 1)", 2)
        push!(w, Xk)
        append!(lbw, [-Inf, -Inf]); append!(ubw, [Inf, Inf])
        push!(G, X_prev - Xk)
    end

    # Allocate an NLP solver
    nlp = Dict("x" => ca.vertcat(w), "f" => J, "g" => ca.vertcat(G))
    solver = ca.nlpsol("solver", "ipopt", nlp)

    # Pass bounds, initial guess and solve NLP
    sol = solver(lbx = lbw, ubx = ubw, lbg = 0.0, ubg = 0.0, x0 = 0.0)

    println("-----")
    println("objective at solution = ", sol["f"])
    println("x0 trajectory = ", sol["x"][1:3:end])
    println("x1 trajectory = ", sol["x"][2:3:end])
end

main()
