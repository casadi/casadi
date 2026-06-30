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
# Demonstration of how the algorithm of an SX function can be accessed and its
# operations traversed.
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

# Create a function
a = ca.SX.sym("a")
b = ca.SX.sym("b", 2)
f = ca.Function("f", [a, b], [2 * a + b], ["a", "b"], ["r"])

# Input values of the same dimensions as the above
input_val = [[2.0], [3.0, 4.0]]

# Output values to be calculated of the same dimensions as the above
output_val = [[0.0, 0.0]]

# Work vector (CasADi addresses it 0-based; we index Julia arrays at +1)
work = zeros(ca.sz_w(f))

function main()
    # Loop over the algorithm
    for k in 0:ca.n_instructions(f) - 1
        op = ca.instruction_id(f, k)        # the atomic operation
        o = ca.instruction_output(f, k)     # output work-vector slots (0-based)
        i = ca.instruction_input(f, k)      # input work-vector slots (0-based)

        if op == ca.OP_CONST
            work[o[1] + 1] = ca.instruction_constant(f, k)
            println("work[ ", o[1], " ] = ", ca.instruction_constant(f, k))
        elseif op == ca.OP_INPUT
            work[o[1] + 1] = input_val[i[1] + 1][i[2] + 1]
            println("work[ ", o[1], " ] = input[ ", i[1], " ][ ", i[2],
                    " ]            ---> ", work[o[1] + 1])
        elseif op == ca.OP_OUTPUT
            output_val[o[1] + 1][o[2] + 1] = work[i[1] + 1]
            println("output[ ", o[1], " ][ ", o[2], " ] = work[ ", i[1],
                    " ]             ---> ", output_val[o[1] + 1][o[2] + 1])
        elseif op == ca.OP_TWICE
            work[o[1] + 1] = 2 * work[i[1] + 1]
            println("work[ ", o[1], " ] = 2*work[ ", i[1],
                    " ]                         ---> ", work[o[1] + 1])
        elseif op == ca.OP_ADD
            work[o[1] + 1] = work[i[1] + 1] + work[i[2] + 1]
            println("work[ ", o[1], " ] = work[ ", i[1], " ] + work[ ", i[2],
                    " ]        ---> ", work[o[1] + 1])
        elseif op == ca.OP_MUL
            work[o[1] + 1] = work[i[1] + 1] * work[i[2] + 1]
            println("work[ ", o[1], " ] = work[ ", i[1], " ] * work[ ", i[2],
                    " ]        ---> ", work[o[1] + 1])
        else
            error("Unknown operation: $op at instruction $k")
        end
    end

    println("------")
    println("Evaluated ", f)
    println("Expected: ", f(ca.DM(input_val[1]), ca.DM(input_val[2])))
    println("Got:      ", output_val[1])
end

main()
