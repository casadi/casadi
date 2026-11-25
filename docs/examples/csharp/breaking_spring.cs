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
using System;
using System.Collections.Generic;
using System.Linq;
using CasADi;

// Simulating and optimizing a spring that breaks with DaeBuilder and event handling
// Joel Andersson, 2024-2025

public static class breaking_spring {
  public static void Main() {
    if (!ca.has_integrator("cvodes")) { Console.WriteLine("SKIP: integrator cvodes not available"); return; }

    // Start with an empty DaeBuilder instance
    var dae = ca.DaeBuilder("breaking_spring");

    // Model variables
    MX m = dae.add("m", "parameter", "tunable", new Dictionary<string, object> { ["start"] = 1.0, ["description"] = "Mass" });
    MX v = dae.add("v", new Dictionary<string, object> { ["start"] = -5.0, ["description"] = "Velocity" });
    MX x = dae.add("x", new Dictionary<string, object> { ["start"] = -1.0, ["description"] = "Displacement" });
    MX k = dae.add("k", "parameter", "tunable", new Dictionary<string, object> { ["start"] = 2.0, ["description"] = "Spring constant" });
    MX c = dae.add("c", "parameter", "tunable", new Dictionary<string, object> { ["start"] = 0.1, ["description"] = "Damping constant" });
    MX d = dae.add("d", "input", new Dictionary<string, object> { ["start"] = 0.0, ["description"] = "Disturbance" });
    MX f = dae.add("f", new Dictionary<string, object> { ["description"] = "Spring force" });
    MX b = dae.add("b", new Dictionary<string, object> { ["type"] = "Boolean", ["start"] = false, ["description"] = "Is the spring broken?" });

    // Equations
    dae.eq(dae.der(x), v);
    dae.eq(f, ca.if_else(ca.logic_not(b), -k * x + d, 0));
    dae.eq(dae.der(v), (f - c * v) / m);
    dae.when(x > 2, [dae.assign("b", 1)]);
    dae.disp(true);

    // Create an integrator in CasADi for hybrid simulation
    double[] tgrid = Enumerable.Range(0, 101).Select(i => 4.0 * i / 100).ToArray();
    var simopts = new Dictionary<string, object> { ["transition"] = dae.transition(), ["verbose"] = false, ["event_tol"] = 1e-12 };
    Function sim = ca.integrator("sim", "cvodes", dae.create(), 0, tgrid, simopts);

    // Disturbance trajectory
    double[] d_traj = new double[tgrid.Length];
    for (int i = 20; i < 30; i++) d_traj[i] = -5;
    for (int i = 40; i < 50; i++) d_traj[i] = 15;

    // Simulate the system
    var res = sim.call(new Dictionary<string, DM> {
      ["x0"] = dae.start(dae.x()).ToArray(), ["p"] = dae.start(dae.p()).ToArray(), ["u"] = d_traj });

    // Visualize the solution
    DM xf = res["xf"].T();
    // plot: states xf (labels dae.x()) over tgrid, and step plot of the disturbance d_traj (label dae.u()[0])
  }

}
