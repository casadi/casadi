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
# Evaluate a function over many inputs, comparing a hand-written vertcat of N
# scalar calls against Function.map(N).
#
# Note: the Python example also tries f.map(N, "openmp"); we use the default
# (serial) map and note it. The expensive function uses 2000 nested sin()
# here (vs 100000) to keep the run snappy.
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

function main()
    # Number of inputs to evaluate
    N = 300

    # Dummy input
    dummyInput = range(0.0, 2.0 * pi, length = N)

    # A moderately expensive function: nested sines
    println("creating dummy function...")
    x = ca.SX.sym("x")
    y = x
    for _ in 1:2000
        y = sin(y)
    end
    f0 = ca.Function("f", [x], [y])

    # Evaluate serially, the old-fashioned way (vertcat of N scalar calls)
    X = ca.MX.sym("x", N)
    Y = ca.vertcat([f0(X[k]) for k in 1:N])
    fNaiveParallel = ca.Function("fParallel", [X], [Y])

    println("evaluating naive parallel function...")
    t0 = time()
    outNaive = fNaiveParallel(ca.DM(collect(dummyInput)))
    println("evaluated naive parallel function in ", round(time() - t0, digits = 3), " seconds")

    # Evaluate it using the serial map construct (row-shaped input/output)
    fMap = f0.map(N)
    println("evaluating serial map function...")
    t0 = time()
    outMap = fMap(ca.DM(collect(dummyInput)'))
    println("evaluated serial map function in ", round(time() - t0, digits = 3), " seconds")

    # The two have differently-shaped outputs (column vs row); compare values
    maxdiff = ca.norm_inf(ca.vec(outNaive) - ca.vec(outMap))
    println("-----")
    println("max |naive - map| = ", maxdiff)
    println("first 5 outputs = ", ca.vec(outMap)[1:5])
end

main()
