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
# Minimal-effort rocket ascent: a single-shooting NLP over the controls of a
# 1-D rocket (position, speed, mass) integrated with explicit Euler.
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

u = ca.MX.sym("u")               # control
x = ca.MX.sym("x", 3)            # state
s = x[1]; v = x[2]; m = x[3]   # position, speed, mass

# ODE right-hand side
xdot = ca.vertcat(v, (u - 0.05 * v * v) / m, -0.1 * u * u)
f = ca.Function("f", [x, u], [xdot])

function main()
    # Integrate with explicit Euler over 0.2 s
    dt = 0.01
    xj = x
    for j in 1:20
        xj = xj + dt * f(xj, u)
    end
    F = ca.Function("F", [x, u], [xj])

    nu = 50                        # number of control segments
    U = ca.MX.sym("U", nu)
    X = ca.MX(ca.DM([0.0, 0.0, 1.0]))
    for k in 1:nu
        X = F(X, U[k])
    end

    J = ca.mtimes(U', U)   # u'*u
    G = X[1:2]
    nlp = Dict("x" => U, "f" => J, "g" => G)

    solver = ca.nlpsol("solver", "ipopt", nlp, Dict("ipopt.tol" => 1e-10, "expand" => true))
    sol = solver(lbx = -0.5, ubx = 0.5, x0 = 0.4, lbg = [10, 0], ubg = [10, 0])

    println("objective at solution = ", sol["f"])
    println("first 5 controls = ", sol["x"][1:5])
end

main()
