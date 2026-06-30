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
# Solve the Rosenbrock problem, formulated as the NLP:
#   minimize     x^2 + 100*z^2
#   subject to   z + (1-x)^2 - y == 0

const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

# Declare variables
x = ca.SX.sym("x")
y = ca.SX.sym("y")
z = ca.SX.sym("z")

# Formulate the NLP
f = x^2 + 100 * z^2
g = z + (1 - x)^2 - y
nlp = Dict("x" => ca.vertcat(x, y, z), "f" => f, "g" => g)

# Create an NLP solver
S = ca.nlpsol("S", "ipopt", nlp)

# Solve the Rosenbrock problem
sol = S(x0 = [2.5, 3.0, 0.75], lbg = 0, ubg = 0)

# Print solution
println()
println("Optimal cost:    ", sol["f"])
println("Primal solution: ", sol["x"])
