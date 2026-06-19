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
# Solve a tiny NLP, then regenerate its solver from self-contained C code:
# `generate_dependencies` emits the NLP functions as C, gcc compiles them into
# a shared library, and a fresh ipopt instance is built straight from that .so.
#
# Unlike the JavaScript port (which only dumps the C text -- the wasm runtime
# has no compiler), the native Julia binding runs the full loop and checks the
# "external" solver reproduces the in-memory solver's solution.
#
#    min x^2 + y^2   s.t.   x + y - 10 = 0
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

function main()
    # Optimization variables
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")

    # Objective and constraint
    f = x * x + y * y
    g = x + y - 10

    # NLP problem structure
    nlp = Dict("x" => ca.vertcat(x, y), "f" => f, "g" => g)

    # Common solver arguments
    arg = (lbx = -Inf, ubx = Inf, lbg = 0, ubg = 0, x0 = 0)

    for mode in ["jit", "external"]
        if mode == "jit"
            # JIT-compile the NLP functions through the "shell" compiler (gcc)
            jit_options = Dict{String,Any}("flags" => ["-O3"], "compiler" => "gcc")
            solver = ca.nlpsol("solver", "ipopt", nlp, Dict("jit" => true, "compiler" => "shell",
                            "jit_options" => jit_options))

        else # mode == "external"
            # Build an ordinary solver, then dump its functions to C and compile
            solver = ca.nlpsol("solver", "ipopt", nlp)
            ca.generate_dependencies(solver, "nlp.c")
            run(`gcc -fPIC -shared -O3 nlp.c -o nlp.so`)

            # Build a new solver instance straight from the compiled code
            solver = ca.nlpsol("solver", "ipopt", "./nlp.so")
        end

        # Solve the NLP
        res = solver(; arg...)

        println("-----")
        println("mode = ", mode)
        println("objective at solution = ", res["f"])
        println("primal solution = ", res["x"])
        println("dual solution (x) = ", res["lam_x"])
        println("dual solution (g) = ", res["lam_g"])
    end
end

main()
