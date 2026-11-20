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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.SX;
import org.casadi.ca;

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
public class dae_reduced_index {
  @SuppressWarnings("unchecked")
  public static void main(String[] args) {
    if (!ca.has_integrator("idas")) { System.out.println("SKIP: integrator idas not available"); return; }
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }

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

    SX alg = ca.vertcat(List.of(dx.minus(u), dy.minus(v), du.minus(T.times(x)), dv.minus(T.times(y)).plus(g),
        x.power(2).plus(y.power(2)).minus(Math.pow(L, 2))));

    // Implicit differential states
    SX x_impl = ca.vertcat(List.of(x, y, u, v));
    // Derivatives of implict differential states
    SX dx_impl = ca.vertcat(List.of(dx, dy, du, dv));

    SX z = T;

    Map<String, SX> dae = Map.of("x_impl", x_impl, "dx_impl", dx_impl, "z", z, "alg", alg);

    // Perform index reduction
    List<Object> r = ca.dae_reduce_index(dae);
    Map<String, SX> dae_reduced = (Map<String, SX>) r.get(0);
    Object stats = r.get(1);
    System.out.println("We had an index 3 system (reduced to index 1 now)");
    System.out.println(stats);
    System.out.println("Here's the reduced DAE:");
    System.out.println(dae_reduced);
    System.out.println("Notice how the algebraic equation now contains a second derivative of x**2+y**2-L**2");
    System.out.println("Also note how x**2+y**2-L**2 and its first derivate are absent from the algebraic equations.");
    System.out.println("Both are kept in 'I' though: the resulting DAE invariants.");

    // The DAE is not yet in a form that CasADi integrator can handle (semi-explicit).
    // Let's convert it, ad obtain some adaptors
    List<Object> m = ca.dae_map_semi_expl(dae, dae_reduced);
    Map<String, SX> dae_se = (Map<String, SX>) m.get(0);
    Function state_to_orig = (Function) m.get(1);
    Function phi = (Function) m.get(2);

    double[] grid = ca.linspace(0, 5, 1000).nonzeros();
    DM tf = new DM(grid).T();
    Function intg = ca.integrator("intg", "idas", dae_se, 0, grid);

    // Before we can integrate, we need a consistant initial guess
    // Let's say we start at y=-0.5, dx=-0.1;
    // consistency equations should automatically identify everything else

    // Encode the desire to only enforce y and dx
    Map<String, DM> init_strength = new LinkedHashMap<>();
    double free = 0; // default
    double force = -1;
    init_strength.put("x_impl", new DM(new double[] {free, force, free, free}));
    init_strength.put("dx_impl", new DM(new double[] {force, free, free, free}));
    init_strength.put("z", new DM(new double[] {free}));

    // Obtain a generating Function for initial guesses
    Function init_gen = ca.dae_init_gen(dae, dae_reduced, "ipopt", init_strength);

    Map<String, Object> init = new LinkedHashMap<>();
    double blank = 0.0; // Will not be enforced but still used as initial guess
    // suggest to initialize with the left-hanging solution by having x=-1 as initial guess
    init.put("x_impl", new double[] {-1.0, -0.5, blank, blank});
    init.put("dx_impl", new double[] {-0.1, blank, blank, blank});
    init.put("z", blank);

    Map<String, GenericMatrix> xz0 = init_gen.invoke(init);
    System.out.println("We have a consistent initial guess now to feed to the integrator:");
    System.out.println(xz0);

    System.out.println("Look, we found that force in the pendulum at t=0 equals:");
    System.out.println(state_to_orig.invoke(Map.of("xf", xz0.get("x0"), "zf", xz0.get("z0"))).get("z"));

    System.out.println("A consistent initial guess should make the invariants zero (up to solver precision):");
    System.out.println(phi.invoke(Map.of("x", xz0.get("x0"), "z", xz0.get("z0"))));

    // Integrate and get resultant xf,zf on grid
    Map<String, GenericMatrix> sol = intg.invoke(xz0);

    DM xf = (DM) sol.get("xf");
    System.out.println("(" + xf.size1() + ", " + xf.size2() + ")");

    // Solution projected into original DAE space
    Map<String, GenericMatrix> sol_orig = state_to_orig.invoke(Map.of("xf", sol.get("xf"), "zf", sol.get("zf")));

    // We can see the large-angle pendulum motion play out well in the u state
    // plot: sol_orig["x_impl"] versus tf, 'Boundary value problem solution trajectory'

    // A perfect integrator will perfectly preserve the values of the invariants over time
    // Integrator errors make the invariants drift in practice
    // This is not a detail; the pendulum length is in fact growing!
    DM error = (DM) phi.invoke(Map.of("x", sol.get("xf"), "z", sol.get("zf"))).get("I");

    // plot: error versus tf, 'Evolution trajectory of invariants'

    // There are some techniques to avoid the drifting of invariants associated with index reduction
    //  - Method of Dummy Derviatives (not implemented)
    //    hard implementation details: may need to switch between different choices online
    //     -> integration in parts triggered by zero-crossing events
    //  - Stop the integration every once in a while and project back (not implemented)
    //  - Baumgarte stabilization (implemented): build into the equations a form of continuous feedback that brings back deviations in invariants back into the origin


    // Demonstrate Baumgarte stabilization with a pole of -1.
    // Drift is absent now.
    r = ca.dae_reduce_index(dae, Map.of("baumgarte_pole", -1));
    dae_reduced = (Map<String, SX>) r.get(0);
    stats = r.get(1);
    m = ca.dae_map_semi_expl(dae, dae_reduced);
    dae_se = (Map<String, SX>) m.get(0);
    state_to_orig = (Function) m.get(1);
    phi = (Function) m.get(2);
    intg = ca.integrator("intg", "idas", dae_se, 0, grid);
    init_gen = ca.dae_init_gen(dae, dae_reduced, "ipopt", init_strength);
    xz0 = init_gen.invoke(init);
    sol = intg.invoke(xz0);
    error = (DM) phi.invoke(Map.of("x", sol.get("xf"), "z", sol.get("zf"))).get("I");
    // plot: error versus tf, 'Boundary value problem solution trajectory with Baumgarte pole=-1'



    // Sweep across different integrator precisions and Baumgarte poles.
    // The choice of pole is not really straightforward

    double[] poles = {-0.1, -1, -10};
    double[] precisions = {0.1, 1, 10};

    double abstol_default = 1e-8;
    double reltol_default = 1e-6;

    for (int i = 0; i < poles.length; i++) {
      double pole = poles[i];
      for (int j = 0; j < precisions.length; j++) {
        double precision = precisions[j];
        r = ca.dae_reduce_index(dae, Map.of("baumgarte_pole", pole));
        dae_reduced = (Map<String, SX>) r.get(0);
        stats = r.get(1);
        m = ca.dae_map_semi_expl(dae, dae_reduced);
        dae_se = (Map<String, SX>) m.get(0);
        state_to_orig = (Function) m.get(1);
        phi = (Function) m.get(2);
        intg = ca.integrator("intg", "idas", dae_se, 0, grid,
            Map.of("abstol", abstol_default * precision, "reltol", reltol_default * precision));
        init_gen = ca.dae_init_gen(dae, dae_reduced, "ipopt", init_strength, Map.of("ipopt.print_level", 0, "print_time", false));
        xz0 = init_gen.invoke(init);
        sol = intg.invoke(xz0);
        long nfevals = ((Number) intg.stats().get("nfevals")).longValue();
        error = (DM) phi.invoke(Map.of("x", sol.get("xf"), "z", sol.get("zf"))).get("I");
        // plot: error versus tf, String.format("pole=%0.1f, abstol=%0.0e => nfevals=%d", pole, abstol_default*precision, nfevals)
      }
    }
  }
}
