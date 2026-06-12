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
# Forward, adjoint and finite-difference sensitivity analysis of an ODE
# integrated over a fixed interval. As in the JS port, the sundials cvodes/idas
# variants are dropped and only rk/collocation are exercised.
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

println("Testing sensitivity analysis in CasADi")
println("******")
println("Testing ODE example")

t = ca.SX.sym("t")               # time
u = ca.SX.sym("u")               # parameter
s = ca.SX.sym("s"); v = ca.SX.sym("v"); m = ca.SX.sym("m")   # differential states
x = ca.vertcat(s, v, m)

alpha = 0.05                   # friction
beta = 0.1                     # fuel consumption rate

ode = ca.vertcat(v, (u - alpha * v * v) / m, -beta * u * u)
quad = v^3 + ((3 - sin(t)) - u)^2
dae = Dict("t" => t, "x" => x, "p" => u, "ode" => ode, "quad" => quad)

tf = 0.5
x0 = ca.DM([0.0, 0.0, 1.0])
u0 = 0.4
h = 0.001

function main()
    for MyIntegrator in ["rk", "collocation"]
        println("========")
        println("Integrator: ", MyIntegrator)
        println("========")

        I = ca.integrator("I", MyIntegrator, dae, 0.0, tf,
                       Dict("number_of_finite_elements" => 100))

        # Unperturbed solution
        res = I(x0 = x0, p = u0)
        xf = res["xf"]; qf = res["qf"]
        println("Unperturbed solution: xf = ", xf, ", qf = ", qf)

        # Finite difference w.r.t. the parameter
        res = I(x0 = x0, p = u0 + h)
        println("Finite difference: d(xf)/d(p) = ", (res["xf"] - xf) / h,
                ", d(qf)/d(p) = ", (res["qf"] - qf) / h)

        # Forward sensitivities via factory
        I_fwd = ca.factory(I, "I_fwd", ["x0", "z0", "p", "fwd:p"], ["fwd:xf", "fwd:qf"])
        res = I_fwd(x0 = x0, p = u0, fwd_p = 1)
        println("Forward: d(xf)/d(p) = ", res["fwd_xf"], ", d(qf)/d(p) = ", res["fwd_qf"])

        # Adjoint sensitivities via factory
        I_adj = ca.factory(I, "I_adj", ["x0", "z0", "p", "adj:qf"], ["adj:x0", "adj:p"])
        res = I_adj(x0 = x0, p = u0, adj_qf = 1)
        adj_x0 = res["adj_x0"]; adj_p = res["adj_p"]
        println("Adjoint: d(qf)/d(x0) = ", adj_x0, ", d(qf)/d(p) = ", adj_p)
    end
end

main()
