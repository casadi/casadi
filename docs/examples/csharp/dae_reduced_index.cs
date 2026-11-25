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
//
// -*- coding: utf-8 -*-
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using CasADi;

// Pendulum: point mass on massless rod
// x'' = T*x, y'' = T*y-g, x**2+y**2-L**2
//
//
//     .
//  L /
//   /
//  /
// o (x,y)
//
public static class dae_reduced_index {
  // Print a dictionary like {key=value, ...}
  static string str(IDictionary d) =>
    "{" + string.Join(", ", d.Keys.Cast<object>().Select(k => $"{k}={d[k]}")) + "}";

  public static void Main() {
    if (!ca.has_integrator("idas")) { Console.WriteLine("SKIP: integrator idas not available"); return; }
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }

    SX x = SX.sym("x"); // horizontal position of point
    SX y = SX.sym("y"); // vertical position of point

    SX dx = SX.sym("dx");
    SX dy = SX.sym("dy");

    SX u = SX.sym("u"); // helper states to bring second order system to first order
    SX v = SX.sym("v");

    SX du = SX.sym("du");
    SX dv = SX.sym("dv");

    // Force
    SX T = SX.sym("T");


    double L = 1;
    double g = 9.81;

    SX alg = ca.vertcat(dx - u, dy - v, du - T * x, dv - T * y + g, ca.pow(x, 2) + ca.pow(y, 2) - Math.Pow(L, 2));

    // Implicit differential states
    SX x_impl = ca.vertcat(x, y, u, v);
    // Derivatives of implict differential states
    SX dx_impl = ca.vertcat(dx, dy, du, dv);

    SX z = T;

    var dae = new Dictionary<string, SX> { ["x_impl"] = x_impl, ["dx_impl"] = dx_impl, ["z"] = z, ["alg"] = alg };

    // Perform index reduction
    var (dae_reduced, stats) = ca.dae_reduce_index(dae);
    Console.WriteLine("We had an index 3 system (reduced to index 1 now)");
    Console.WriteLine(str(stats));
    Console.WriteLine("Here's the reduced DAE:");
    Console.WriteLine(str(dae_reduced));
    Console.WriteLine("Notice how the algebraic equation now contains a second derivative of x**2+y**2-L**2");
    Console.WriteLine("Also note how x**2+y**2-L**2 and its first derivate are absent from the algebraic equations.");
    Console.WriteLine("Both are kept in 'I' though: the resulting DAE invariants.");

    // The DAE is not yet in a form that CasADi integrator can handle (semi-explicit).
    // Let's convert it, ad obtain some adaptors
    var (dae_se, state_to_orig, phi) = ca.dae_map_semi_expl(dae, dae_reduced);

    List<double> grid = ca.linspace(0.0, 5.0, 1000).nonzeros();
    DM tf = ((DM)grid.ToArray()).T();
    Function intg = ca.integrator("intg", "idas", dae_se, 0, grid);

    // Before we can integrate, we need a consistant initial guess
    // Let's say we start at y=-0.5, dx=-0.1;
    // consistency equations should automatically identify everything else

    // Encode the desire to only enforce y and dx
    double free = 0; // default
    double force = -1;
    var init_strength = new Dictionary<string, DM> {
      ["x_impl"] = new double[] { free, force, free, free },
      ["dx_impl"] = new double[] { force, free, free, free },
      ["z"] = new double[] { free } };

    // Obtain a generating Function for initial guesses
    Function init_gen = ca.dae_init_gen(dae, dae_reduced, "ipopt", init_strength);

    double blank = 0.0; // Will not be enforced but still used as initial guess
    // suggest to initialize with the left-hanging solution by having x=-1 as initial guess
    var init = new Dictionary<string, DM> {
      ["x_impl"] = new double[] { -1.0, -0.5, blank, blank },
      ["dx_impl"] = new double[] { -0.1, blank, blank, blank },
      ["z"] = blank };

    var xz0 = init_gen.call(init);
    Console.WriteLine("We have a consistent initial guess now to feed to the integrator:");
    Console.WriteLine(str(xz0));

    Console.WriteLine("Look, we found that force in the pendulum at t=0 equals:");
    Console.WriteLine(state_to_orig.call(new Dictionary<string, DM> { ["xf"] = xz0["x0"], ["zf"] = xz0["z0"] })["z"]);

    Console.WriteLine("A consistent initial guess should make the invariants zero (up to solver precision):");
    Console.WriteLine(str(phi.call(new Dictionary<string, DM> { ["x"] = xz0["x0"], ["z"] = xz0["z0"] })));

    // Integrate and get resultant xf,zf on grid
    var sol = intg.call(xz0);

    DM xf = sol["xf"];
    Console.WriteLine($"({xf.size1()}, {xf.size2()})");

    // Solution projected into original DAE space
    var sol_orig = state_to_orig.call(new Dictionary<string, DM> { ["xf"] = sol["xf"], ["zf"] = sol["zf"] });

    // We can see the large-angle pendulum motion play out well in the u state
    // plot: sol_orig["x_impl"] versus tf, 'Boundary value problem solution trajectory'

    // A perfect integrator will perfectly preserve the values of the invariants over time
    // Integrator errors make the invariants drift in practice
    // This is not a detail; the pendulum length is in fact growing!
    DM error = phi.call(new Dictionary<string, DM> { ["x"] = sol["xf"], ["z"] = sol["zf"] })["I"];

    // plot: error versus tf, 'Evolution trajectory of invariants'

    // There are some techniques to avoid the drifting of invariants associated with index reduction
    //  - Method of Dummy Derviatives (not implemented)
    //    hard implementation details: may need to switch between different choices online
    //     -> integration in parts triggered by zero-crossing events
    //  - Stop the integration every once in a while and project back (not implemented)
    //  - Baumgarte stabilization (implemented): build into the equations a form of continuous feedback that brings back deviations in invariants back into the origin


    // Demonstrate Baumgarte stabilization with a pole of -1.
    // Drift is absent now.
    (dae_reduced, stats) = ca.dae_reduce_index(dae, new Dictionary<string, object> { ["baumgarte_pole"] = -1.0 });
    (dae_se, state_to_orig, phi) = ca.dae_map_semi_expl(dae, dae_reduced);
    intg = ca.integrator("intg", "idas", dae_se, 0, grid);
    init_gen = ca.dae_init_gen(dae, dae_reduced, "ipopt", init_strength);
    xz0 = init_gen.call(init);
    sol = intg.call(xz0);
    error = phi.call(new Dictionary<string, DM> { ["x"] = sol["xf"], ["z"] = sol["zf"] })["I"];
    // plot: error versus tf, 'Boundary value problem solution trajectory with Baumgarte pole=-1'



    // Sweep across different integrator precisions and Baumgarte poles.
    // The choice of pole is not really straightforward

    double[] poles = { -0.1, -1, -10 };
    double[] precisions = { 0.1, 1, 10 };

    double abstol_default = 1e-8;
    double reltol_default = 1e-6;

    foreach (double pole in poles) {
      foreach (double precision in precisions) {
        (dae_reduced, stats) = ca.dae_reduce_index(dae, new Dictionary<string, object> { ["baumgarte_pole"] = pole });
        (dae_se, state_to_orig, phi) = ca.dae_map_semi_expl(dae, dae_reduced);
        intg = ca.integrator("intg", "idas", dae_se, 0, grid,
            new Dictionary<string, object> { ["abstol"] = abstol_default * precision, ["reltol"] = reltol_default * precision });
        init_gen = ca.dae_init_gen(dae, dae_reduced, "ipopt", init_strength,
            new Dictionary<string, object> { ["ipopt.print_level"] = 0L, ["print_time"] = false });
        xz0 = init_gen.call(init);
        sol = intg.call(xz0);
        long nfevals = Convert.ToInt64(intg.stats()["nfevals"]);
        error = phi.call(new Dictionary<string, DM> { ["x"] = sol["xf"], ["z"] = sol["zf"] })["I"];
        // plot: error versus tf, $"pole={pole:0.0}, abstol={abstol_default * precision:0e0} => nfevals={nfevals}"
      }
    }
  }
}
