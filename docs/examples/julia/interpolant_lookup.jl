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
# Interpolant lookup tables (no Python/JavaScript twin).
#
# CasADi's `interpolant` turns sampled data into a Function that is fully
# differentiable: a data-driven lookup table behaves like any other node in the
# expression graph, so you can embed it in an NLP/integrator and let AD flow
# through it. This example builds a 1-D bspline LUT and a 2-D linear LUT, then
# differentiates both symbolically.
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

function main()
    # ---- 1-D smooth lookup table: a bspline through samples of f(x) = x^2 ----
    xgrid = collect(range(0.0, 6.0, length = 7))     # knots [0, 1, ..., 6]
    yvals = xgrid .^ 2                               # the sampled values
    LUT = ca.interpolant("LUT", "bspline", [xgrid], yvals)

    println("1-D bspline LUT (samples of x^2)")
    println("  LUT(2.5) = ", LUT(2.5), "   (exact 2.5^2 = ", 2.5^2, ")")
    println("  LUT(4.0) = ", LUT(4.0))

    # Embed the table in a symbolic expression and differentiate through it.
    x = ca.MX.sym("x")
    expr = LUT(x) + sin(x)
    dexpr = ca.Function("dexpr", [x], [ca.gradient(expr, x)])
    println("  d/dx[LUT(x) + sin(x)] at x = 3 : ", dexpr(3.0),
            "   (exact 2*3 + cos(3) = ", 2 * 3 + cos(3.0), ")")

    # ---- 2-D linear lookup table on a non-uniform grid ----
    grid = [[0.0, 1, 4, 5], [0.0, 2, 3]]            # x-knots, y-knots
    values = [0.0, 1, 8, 3,                         # row-major over the grid
              10, -11, 12, 13,
              20, 31, -42, 53]
    F = ca.interpolant("F", "linear", grid, values)

    println("2-D linear LUT on a 4x3 grid")
    println("  F([1, 2])   = ", F(ca.DM([1.0, 2.0])), "   (a grid node, exact -11)")
    println("  F([1, 2.4]) = ", F(ca.DM([1.0, 2.4])),
            "   (interpolated, exact ", -11 + 0.4 * (31 + 11), ")")

    # Jacobian of the 2-D table wrt its query point (bilinear -> piecewise const)
    X = ca.MX.sym("X", 2)
    J = ca.Function("J", [X], [ca.jacobian(F(X), X)])
    println("  dF/dX at [3, 2.4] = ", J(ca.DM([3.0, 2.4])))
end

main()
