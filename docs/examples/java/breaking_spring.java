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
import org.casadi.*;
import java.util.List;
import java.util.Map;

// Simulating and optimizing a spring that breaks with DaeBuilder and event handling
// Joel Andersson, 2024-2025

public class breaking_spring {
  public static void main(String[] args) {
    if (!ca.has_integrator("cvodes")) { System.out.println("SKIP: integrator cvodes not available"); return; }

    // Start with an empty DaeBuilder instance
    DaeBuilder dae = new DaeBuilder("breaking_spring");

    // Model variables
    MX m = dae.add("m", "parameter", "tunable", Map.of("start", 1, "description", "Mass"));
    MX v = dae.add("v", Map.of("start", -5, "description", "Velocity"));
    MX x = dae.add("x", Map.of("start", -1, "description", "Displacement"));
    MX k = dae.add("k", "parameter", "tunable", Map.of("start", 2, "description", "Spring constant"));
    MX c = dae.add("c", "parameter", "tunable", Map.of("start", 0.1, "description", "Damping constant"));
    MX d = dae.add("d", "input", Map.of("start", 0, "description", "Disturbance"));
    MX f = dae.add("f", Map.of("description", "Spring force"));
    MX b = dae.add("b", Map.of("type", "Boolean", "start", false, "description", "Is the spring broken?"));

    // Equations
    dae.eq(dae.der(x), v);
    dae.eq(f, ca.if_else(ca.logic_not(b), k.unaryMinus().times(x).plus(d), 0));
    dae.eq(dae.der(v), f.minus(c.times(v)).div(m));
    dae.when(x.gt(2), List.of(dae.assign("b", 1)));
    dae.disp(true);

    // Create an integrator in CasADi for hybrid simulation
    double[] tgrid = new double[101];
    for (int i = 0; i < 101; i++) tgrid[i] = 4.0 * i / 100;
    Map<String, Object> simopts = Map.of("transition", dae.transition(), "verbose", false, "event_tol", 1e-12);
    Function sim = ca.integrator("sim", "cvodes", dae.create(), 0, tgrid, simopts);

    // Disturbance trajectory
    double[] d_traj = new double[tgrid.length];
    for (int i = 20; i < 30; i++) d_traj[i] = -5;
    for (int i = 40; i < 50; i++) d_traj[i] = 15;

    // Simulate the system
    Map<String, GenericMatrix> res = sim.invoke(Map.of("x0", dae.start(dae.x()), "p", dae.start(dae.p()), "u", d_traj));

    // Visualize the solution
    DM xf = ((DM) res.get("xf")).T();
    // plot: states xf (labels dae.x()) over tgrid, and step plot of the disturbance d_traj (label dae.u().get(0))
  }

}
