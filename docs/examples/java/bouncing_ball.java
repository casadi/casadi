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

public class bouncing_ball {
  public static void main(String[] args) {
    if (!ca.has_integrator("cvodes")) { System.out.println("SKIP: integrator cvodes not available"); return; }

    // Height and velocity of ball
    SX h = SX.sym("h");
    SX v = SX.sym("v");
    SX x = ca.vertcat(List.of(h, v));

    // ODE right-hand-side
    SX hdot = v;
    SX vdot = new SX(-9.81);
    SX xdot = ca.vertcat(List.of(hdot, vdot));

    // Event indicator, trigger when it becomes negative
    SX event_indicator = h;

    // DAE problem structure, with zero-crossing output
    Map<String, SX> dae = Map.of("x", x, "ode", xdot, "zero", event_indicator);

    // Event transition function
    SX post_x = ca.vertcat(List.of(h, v.times(-0.8)));
    Function transition = new Function("transition", Map.of("x", x, "post_x", post_x),
                                       ca.event_in(), ca.event_out());

    // Create an integrator instance for integrating over 7s
    double[] tgrid = new double[100];
    for (int i = 0; i < 100; i++) tgrid[i] = 7.0 * i / 99;
    Function sim = ca.integrator("sim", "cvodes", dae, 0, tgrid,
                                 Map.of("transition", transition));

    // Simulate with initial height of 5
    double[] x0 = {5, 0};
    Map<String, GenericMatrix> simres = sim.invoke(Map.of("x0", x0));

    // Visualize the solution
    DM h_traj = ((DM) simres.get("xf")).get(0, new Slice()).T();
    // plot: h_traj over tgrid
  }
}
