//
//     MIT No Attribution
//
//     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//     Permission is hereby granted, free of charge, to any person obtaining a copy of this
//     software and associated documentation files (the "Software"), to deal in the Software
//     without restriction, including without limitation the rights to use, copy, modify,
//     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//     permit persons to whom the Software is furnished to do so.
//
//     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
// Interpolant lookup tables (no Python/JavaScript twin).
//
// CasADi's `interpolant` turns sampled data into a Function that is fully
// differentiable: a data-driven lookup table behaves like any other node in the
// expression graph, so you can embed it in an NLP/integrator and let AD flow
// through it. This example builds a 1-D bspline LUT and a 2-D linear LUT, then
// differentiates both symbolically.
using System;
using System.Linq;
using CasADi;

public static class interpolant_lookup {
  public static void Main() {
    if (!ca.has_interpolant("bspline")) { Console.WriteLine("SKIP: interpolant bspline not available"); return; }
    if (!ca.has_interpolant("linear")) { Console.WriteLine("SKIP: interpolant linear not available"); return; }

    // ---- 1-D smooth lookup table: a bspline through samples of f(x) = x^2 ----
    double[] xgrid = Enumerable.Range(0, 7).Select(i => (double)i).ToArray(); // knots [0, 1, ..., 6]
    double[] yvals = xgrid.Select(v => v * v).ToArray();                       // the sampled values
    Function LUT = ca.interpolant("LUT", "bspline", [xgrid], yvals);

    Console.WriteLine("1-D bspline LUT (samples of x^2)");
    Console.WriteLine($"  LUT(2.5) = {LUT.Call(2.5)}   (exact 2.5^2 = {Math.Pow(2.5, 2)})");
    Console.WriteLine($"  LUT(4.0) = {LUT.Call(4.0)}");

    // Embed the table in a symbolic expression and differentiate through it.
    MX x = MX.sym("x");
    MX expr = LUT.Call(x) + ca.sin(x);
    var dexpr = ca.Function("dexpr", [x], [ca.gradient(expr, x)]);
    Console.WriteLine($"  d/dx[LUT(x) + sin(x)] at x = 3 : {dexpr.Call(3.0)}"
        + $"   (exact 2*3 + cos(3) = {2 * 3 + Math.Cos(3.0)})");

    // ---- 2-D linear lookup table on a non-uniform grid ----
    double[][] grid = { new[] { 0.0, 1, 4, 5 }, [0.0, 2, 3] };  // x-knots, y-knots
    double[] values = { 0.0, 1, 8, 3,                                   // row-major over the grid
                        10, -11, 12, 13,
                        20, 31, -42, 53 };
    Function F = ca.interpolant("F", "linear", grid, values);

    Console.WriteLine("2-D linear LUT on a 4x3 grid");
    Console.WriteLine($"  F([1, 2])   = {F.Call(new double[] { 1.0, 2.0 })}   (a grid node, exact -11)");
    Console.WriteLine($"  F([1, 2.4]) = {F.Call(new double[] { 1.0, 2.4 })}"
        + $"   (interpolated, exact {-11 + 0.4 * (31 + 11)})");

    // Jacobian of the 2-D table wrt its query point (bilinear -> piecewise const)
    MX X = MX.sym("X", 2);
    var J = ca.Function("J", [X], [ca.jacobian(F.Call(X), X)]);
    Console.WriteLine($"  dF/dX at [3, 2.4] = {J.Call(new double[] { 3.0, 2.4 })}");
  }
}
