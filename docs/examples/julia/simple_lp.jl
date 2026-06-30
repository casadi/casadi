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
# A linear program is a Conic (QP) problem with no quadratic term: pass only
# the sparsity of the constraint matrix `a`.
#
#   minimize    3x + 4y
#   subject to  x + 2y <= 14
#               3x -  y >= 0
#                x -  y <= 2
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

# Sparsity of the LP linear term (3 constraints x 2 variables)
A = ca.dense(ca.Sparsity, 3, 2)

# Create solver (qpoases handles LPs as a degenerate QP)
solver = ca.conic("solver", "qpoases", Dict("a" => A))

g = ca.DM([3.0, 4.0])
a = ca.DM([1.0 2.0; 3.0 -1.0; 1.0 -1.0])
lba = ca.DM([-Inf, 0.0, -Inf])
uba = ca.DM([14.0, Inf, 2.0])

sol = solver(g = g, a = a, lba = lba, uba = uba)
for k in keys(sol)
    println(k, " = ", sol[k])
end
