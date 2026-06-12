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
# Generate C source for the gradient of a 7x7 determinant, compile it with gcc
# at several optimization levels, load each shared library back through
# `external`, and check the generated code reproduces the original numerically.
#
# Unlike the JavaScript port (which only dumps the C text -- the wasm runtime
# has no compiler), the native Julia binding runs the full codegen loop:
# generate -> gcc -> external -> call.
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

function main()
    # Form an expression for the gradient of the determinant
    x = ca.SX.sym("x", 7, 7)
    gd = ca.gradient(ca.det(x), x)

    # Random point to evaluate it
    x0 = ca.rand(ca.DM, 7, 7)

    # Form a function and generate C code
    name = "grad_det"
    grad_det = ca.Function(name, [x], [gd], ["x"], ["gd"])
    cname = ca.generate(grad_det, name * ".c")
    println("generated C source: ", cname)

    # Compile the generated C at three optimization levels
    flags = Dict("no_opt" => String[], "O3" => ["-O3"], "Os" => ["-Os"])
    sofiles = Dict{String,String}()
    for (label, fl) in flags
        oname = name * "_" * label * ".so"
        print("compiling ", oname, " ... ")
        t1 = time()
        run(`gcc -fPIC -shared $fl $cname -o $oname`)
        println(round((time() - t1) * 1e3, digits = 1), " ms")
        sofiles[label] = oname
    end

    # Read each compiled function back in via external()
    f_test = [("symbolic", grad_det)]
    for label in ("no_opt", "O3", "Os")
        push!(f_test, (label, ca.external(name, "./" * sofiles[label])))
    end

    # The reference result is the symbolic function evaluated at x0
    ref = ca.nonzeros(grad_det(x0))
    num_op = ca.n_nodes(grad_det)

    for (label, f) in f_test
        nrep = 10000
        t1 = time()
        local r
        for _ in 1:nrep
            r = f(x0)
        end
        dt = (time() - t1) / nrep
        rz = ca.nonzeros(r)
        # The compiled code must reproduce the symbolic result bit-for-bit
        @assert maximum(abs.(rz .- ref)) < 1e-12
        println("[", label, "] result[1:3] = ", round.(rz[1:3], digits = 6),
                "  time = ", round(dt * 1e3, digits = 4), " ms")
    end

    println("number of elementary operations: ", num_op)
    println("all compiled variants match the symbolic gradient.")
end

main()
