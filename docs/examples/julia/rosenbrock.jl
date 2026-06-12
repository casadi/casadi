#
#     MIT No Attribution
#
#     Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
#
# Julia port of docs/examples/python/rosenbrock.py.
#
# Solve the Rosenbrock problem, formulated as the NLP:
#   minimize     x^2 + 100*z^2
#   subject to   z + (1-x)^2 - y == 0

const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
using .CasADiNative

# Declare variables
x = sym(SX, "x")
y = sym(SX, "y")
z = sym(SX, "z")

# Formulate the NLP
f = x^2 + 100 * z^2
g = z + (1 - x)^2 - y
nlp = Dict("x" => vertcat(x, y, z), "f" => f, "g" => g)

# Create an NLP solver
S = nlpsol("S", "ipopt", nlp)

# Solve the Rosenbrock problem
sol = S(x0 = [2.5, 3.0, 0.75], lbg = 0, ubg = 0)

# Print solution
println()
println("Optimal cost:    ", sol["f"])
println("Primal solution: ", sol["x"])
