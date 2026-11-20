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
import java.util.List;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.MX;
import org.casadi.ca;

public class interpolant_lookup {
  public static void main(String[] args) {
    if (!ca.has_interpolant("bspline")) { System.out.println("SKIP: interpolant bspline not available"); return; }
    if (!ca.has_interpolant("linear")) { System.out.println("SKIP: interpolant linear not available"); return; }

    // ---- 1-D smooth lookup table: a bspline through samples of f(x) = x^2 ----
    double[] xgrid = new double[7];                   // knots [0, 1, ..., 6]
    for (int i = 0; i < 7; i++) xgrid[i] = i;
    double[] yvals = new double[7];                   // the sampled values
    for (int i = 0; i < 7; i++) yvals[i] = xgrid[i] * xgrid[i];
    Function LUT = ca.interpolant("LUT", "bspline", new double[][] {xgrid}, yvals);

    System.out.println("1-D bspline LUT (samples of x^2)");
    System.out.println("  LUT(2.5) = " + LUT.invoke(2.5) + "   (exact 2.5^2 = " + Math.pow(2.5, 2) + ")");
    System.out.println("  LUT(4.0) = " + LUT.invoke(4.0));

    // Embed the table in a symbolic expression and differentiate through it.
    MX x = MX.sym("x");
    MX expr = ((MX) LUT.invoke(x)).plus(ca.sin(x));
    Function dexpr = new Function("dexpr", List.of(x), List.of(ca.gradient(expr, x)));
    System.out.println("  d/dx[LUT(x) + sin(x)] at x = 3 : " + dexpr.invoke(3.0)
        + "   (exact 2*3 + cos(3) = " + (2 * 3 + Math.cos(3.0)) + ")");

    // ---- 2-D linear lookup table on a non-uniform grid ----
    double[][] grid = {{0.0, 1, 4, 5}, {0.0, 2, 3}};  // x-knots, y-knots
    double[] values = {0.0, 1, 8, 3,                  // row-major over the grid
                       10, -11, 12, 13,
                       20, 31, -42, 53};
    Function F = ca.interpolant("F", "linear", grid, values);

    System.out.println("2-D linear LUT on a 4x3 grid");
    System.out.println("  F([1, 2])   = " + F.invoke(new DM(new double[] {1.0, 2.0})) + "   (a grid node, exact -11)");
    System.out.println("  F([1, 2.4]) = " + F.invoke(new DM(new double[] {1.0, 2.4}))
        + "   (interpolated, exact " + (-11 + 0.4 * (31 + 11)) + ")");

    // Jacobian of the 2-D table wrt its query point (bilinear -> piecewise const)
    MX X = MX.sym("X", 2);
    Function J = new Function("J", List.of(X), List.of(ca.jacobian((MX) F.invoke(X), X)));
    System.out.println("  dF/dX at [3, 2.4] = " + J.invoke(new DM(new double[] {3.0, 2.4})));
  }
}
