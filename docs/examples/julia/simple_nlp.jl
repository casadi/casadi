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
# Solve a small NLP:
#   minimize     x0^2 + x1^2
#   subject to   x0 + x1 - 10 >= 0
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

# Declare variables
x = ca.SX.sym("x", 2)

# Form the NLP
f = x[1]^2 + x[2]^2          # objective
g = x[1] + x[2] - 10         # constraint
nlp = Dict("x" => x, "f" => f, "g" => g)

# Allocate an ipopt solver and solve with g >= 0
solver = ca.nlpsol("solver", "ipopt", nlp)
sol = solver(lbg = 0)

# Print solution
println("-----")
println("objective at solution = ", sol["f"])
println("primal solution = ", sol["x"])
println("dual solution (x) = ", sol["lam_x"])
println("dual solution (g) = ", sol["lam_g"])
