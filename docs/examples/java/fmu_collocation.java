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

// Example on how to formulate and solve an optimal control problem using direct collocation.
// The model is a semi-explicit DAE.
//
// Based on collocation example by Mario Zanon, Sebastien Gross 2012
// Reformulated to a semi-explict formulation and modern syntax by Joel Andersson, 2026

import org.casadi.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.zip.*;

public class fmu_collocation {

  // Python repr of a str -> str dict
  static String pyrepr(Map<String, ?> m) {
    StringJoiner s = new StringJoiner(", ", "{", "}");
    for (Map.Entry<String, ?> e : m.entrySet()) s.add("'" + e.getKey() + "': '" + e.getValue() + "'");
    return s.toString();
  }

  static double[] select(List<String> names, List<String> which, double value) {
    double[] r = new double[names.size()];
    for (int i = 0; i < r.length; i++) r[i] = which.contains(names.get(i)) ? value : 0;
    return r;
  }

  @SuppressWarnings("unchecked")
  public static void main(String[] args) throws IOException {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }
    double inf = Double.POSITIVE_INFINITY;

    // Precompiled FMU?
    String fmu = null; // Generate CasADi FMU
    // fmu = "casadi_symbolic_fmu"; // FMU from CasADi, use serialized CasADi expressions (default)
    // fmu = "casadi_blackbox_fmu"; // FMU from CasADi, do not use serialized CasADi expressions
    // fmu = "dymola_fmu"; // FMU from Dymola, compiled from Modelica

    DaeBuilder dae;
    if (fmu == null) {
      // Start with an empty DaeBuilder instance
      dae = new DaeBuilder("Crane");

      // Hard coded constants
      double mL = 1.; // Mass of the load
      double mC = 1.; // Mass of the crane
      double g = 9.81; // Gravity
      double length = 1.0; // Rope/pole length
      boolean inverted = true; // Inverted pendulum on a cart rather than crane with hanging load

      // States
      MX xL = dae.add("xL", Map.of("initial", "exact", "start", 0,
          "description", "Horizontal position of the load"));
      MX y = dae.add("y", Map.of("initial", "exact", "start", inverted ? length : -length,
          "description", "Vertical position of the load"));
      MX xC = dae.add("xC", Map.of("initial", "exact", "start", 0,
          "description", "Horizontal position of the crane"));
      MX dxL = dae.add("dxL", Map.of("initial", "exact", "start", 0,
          "description", "Horizontal velocity of the load"));
      MX dy = dae.add("dy", Map.of("initial", "exact", "start", 0,
          "description", "Vertical velocity of the load"));
      MX dxC = dae.add("dxC", Map.of("initial", "exact", "start", 0,
          "description", "Horizontal velocity of the crane"));

      // Input
      MX u = dae.add("u", "input", Map.of("min", -2, "max", 2, "start", 0,
          "description", "Horizontal force applied to the crane"));

      // Algebraic variables
      MX T = dae.add("T", Map.of("initial", "approx", "start", inverted ? g : -g,
          "description", "Tension in the rope/pole"));

      // Ordinary differential equations
      MX ddxL = xL.minus(xC).div(length).times(T).div(mL);
      MX ddy = ca.minus(g, y.div(length).times(T).div(mL));
      MX ddxC = xL.minus(xC).div(length).times(T).minus(u).div(mC);
      dae.eq(dae.der(xL), dxL);
      dae.eq(dae.der(y), dy);
      dae.eq(dae.der(xC), dxC);
      dae.eq(dae.der(dxL), ddxL);
      dae.eq(dae.der(dy), ddy);
      dae.eq(dae.der(dxC), ddxC);

      // Original algebraic equation
      MX res = ca.power(y, 2).plus(ca.power(xL.minus(xC), 2)).minus(Math.pow(length, 2));
      System.out.println("Original residual equation: " + res);

      // Index reduction: Differentiate with respect to time
      // In CasADi, the most efficient way to do index reduction it to use forward
      // mode AD, e.g. us using Jacobian-times-vector products (jtimes):
      // Using the chain rule we have: d(res)/dt = d(res)/dx * dx/dt
      MX res1 = ca.jtimes(res, ca.vertcat(dae.inputs("x")), ca.vertcat(dae.outputs("ode")));
      MX res2 = ca.jtimes(res1, ca.vertcat(dae.inputs("x")), ca.vertcat(dae.outputs("ode")));
      System.out.println("Residual equation after index reduction: " + res2);

      // Add algebraic equations
      dae.eq(0, res2);

      // Simulation duration
      dae.set_start_time(0);
      dae.set_stop_time(5);

      // Generate FMU
      Map<String, Object> fmu_files = new LinkedHashMap<>(new TreeMap<>(dae.export_fmu()));
      System.out.println("Generated files: " + pyrepr(fmu_files));

      // Compile DLL
      Path fmi_headers = Paths.get("../../../external_packages/FMI-Standard-3.0/headers");
      StringJoiner cfiles = new StringJoiner(" ");
      for (String f : fmu_files.keySet()) if (f.endsWith(".c")) cfiles.add(f);
      // Compile sources (in Linux)
      if (File.separatorChar == '/') {
        try {
          String sofile = dae.name() + ".so";
          new ProcessBuilder("sh", "-c", "gcc --shared -fPIC -I" + fmi_headers + " " + cfiles + " -o " + sofile)
              .inheritIO().start().waitFor();
          System.out.println("Compiled " + sofile);
          fmu_files.put(sofile, "binaries/x86_64-linux/" + sofile);
        } catch (IOException | InterruptedException e) {
          System.out.println("Failed to create FMU binary: " + e);
        }
      }

      // Package into an FMU
      String fmuname = dae.name() + ".fmu";
      try (ZipOutputStream fmufile = new ZipOutputStream(new FileOutputStream(fmuname))) {
        for (Map.Entry<String, Object> e : fmu_files.entrySet()) {
          fmufile.putNextEntry(new ZipEntry(Paths.get((String) e.getValue()).normalize().toString()));
          Files.copy(Paths.get(e.getKey()), fmufile);
          fmufile.closeEntry();
          Files.delete(Paths.get(e.getKey()));
        }
      }
      System.out.println("Created FMU: " + fmuname);
    } else if (fmu.equals("casadi_symbolic_fmu")) {
      // Use FMU generated from CasADi using the script above: Standard FMU + CasADi expressions + LS-DAE
      dae = new DaeBuilder("Crane", "./Crane.fmu");
    } else if (fmu.equals("casadi_blackbox_fmu")) {
      // Use FMU generated from CasADi using the script above: Standard FMU + LS-DAE
      dae = new DaeBuilder("Crane", "./Crane.fmu", Map.of("enable_ls_serialization", false));
    } else if (fmu.equals("dymola_fmu")) {
      // Use FMU compiled from Modelica using Dymola
      dae = new DaeBuilder("crane", "../../../test/data/Crane.fmu");
    } else {
      throw new IllegalArgumentException("Unknown FMU option: " + fmu);
    }

    // Display the model
    dae.disp(true);

    // DAE right-hand-side function
    Function ffcn = dae.create("ffcn", List.of("x", "z", "u"), List.of("ode", "alg", "y"));

    // Variable names
    List<String> U = dae.u();
    List<String> X = dae.x();
    List<String> Z = dae.z();
    List<String> Y = dae.y();

    // Refence location for crane and load
    double xref = 0.1;

    // Quadratic cost at the beginning of control intervals
    DM X_ref = new DM(select(X, List.of("xL", "xC"), xref));
    DM Y_ref = new DM(select(Y, List.of("xL", "xC"), xref));
    DM X_weights = new DM(select(X, List.of("xL", "xC"), 1.));
    DM Y_weights = new DM(select(Y, List.of("xL", "xC"), 1.));

    // Quadratic cost at final time
    DM Xf_ref = new DM(select(X, List.of("xL", "xC"), xref));
    DM Yf_ref = new DM(select(Y, List.of("xL", "xC"), xref));
    DM Xf_weights = new DM(select(X, List.of("xL", "xC", "dxL", "dy"), 1.));
    DM Yf_weights = new DM(select(Y, List.of("xL", "xC", "dxL", "dy"), 1.));

    // Discretization
    int nk = 50; // Number of control intervals
    int deg = 4; // Degree of interpolating polynomial

    // Size of the finite elements
    double tf = 5.0;
    double h = tf / nk;

    // We will use collocation discretization using Legendre roots, cf.
    // Nonlinear Programming: Concepts, Algorithms, and Applications to Chemical Processes
    // by Lorenz Biegler (2010).
    // The roots can be queried from CasADi or looked up in the above textbook
    double[] tau_root = ca.collocation_points(deg, "radau");

    // We can query the coefficient for the interpolating polynomials using a helper function
    // C: Coefficients for calculating state derivatives
    // D: Coefficients for calculating state at the end of the interval
    // B: Coefficients for calculating a quadrature at the end of the interval
    List<DM> CDB = ca.collocation_coeff(tau_root);
    DM C = CDB.get(0), D = CDB.get(1), B = CDB.get(2);

    // We will construct an NLP incrementally, starting with no decision variables,
    // and constraints, and zero objective
    List<MX> W = new ArrayList<>();
    List<MX> G = new ArrayList<>();
    MX J = new MX(0);

    // Solution trajectories
    List<MX> U_all = new ArrayList<>();
    List<MX> X_all = new ArrayList<>();
    List<MX> Z_all = new ArrayList<>();
    List<MX> Y_all = new ArrayList<>();

    // Time grids, for solution
    List<Double> Ugrid = new ArrayList<>();
    List<Double> Xgrid = new ArrayList<>();
    List<Double> Zgrid = new ArrayList<>();
    List<Double> Ygrid = new ArrayList<>();

    // Current time
    double tk = 0;

    // Add a variable for the initial state
    MX Xk = MX.sym("X0", X.size());
    X_all.add(Xk);
    W.add(Xk);
    Xgrid.add(tk);
    MX Yk = null;
    // Loop over all control intervals, constructing the NLP and x, u trajectories
    for (int k = 0; k < nk; k++) {
      // Add contribution to the cost
      MX Dx = Xk.minus(X_ref);
      J = J.plus(ca.dot(Dx, Dx.times(X_weights)));
      if (k > 0 && Y.size() > 0) {
        // Add contribution to the cost
        MX Dy = Yk.minus(Y_ref);
        J = J.plus(ca.dot(Dy, Dy.times(Y_weights)));
      }
      // Symbolic expression for the control for the interval
      MX Uk = MX.sym("u" + k, U.size());
      U_all.add(Uk);
      W.add(Uk);
      Ugrid.add(tk);
      // Declare variables at each collocation point
      List<MX> Xc = new ArrayList<>();
      List<MX> Zc = new ArrayList<>();
      List<MX> Yc = new ArrayList<>();
      for (int j = 1; j < deg + 1; j++) {
        // Differential state at collocation points
        MX Xkj = MX.sym("X" + k + "_" + j, X.size());
        Xc.add(Xkj);
        X_all.add(Xkj);
        W.add(Xkj);
        Xgrid.add(tk + tau_root[j - 1] * h);
        // Algebraic variables at collocation points
        MX Zkj = MX.sym("Z" + k + "_" + j, Z.size());
        Zc.add(Zkj);
        Z_all.add(Zkj);
        W.add(Zkj);
        Zgrid.add(tk + tau_root[j - 1] * h);
        // Output variables at collocation points
        MX Ykj = MX.sym("Y" + k + "_" + j, Y.size());
        Yc.add(Ykj);
        Y_all.add(Ykj);
        W.add(Ykj);
        Ygrid.add(tk + tau_root[j - 1] * h);
      }
      // Collect the collocation equations
      for (int j = 0; j < deg; j++) {
        // Expression for the state derivative at the collocation point
        MX Xdot_j = Xk.times(C.get(0, j));
        for (int r = 0; r < deg; r++) Xdot_j = Xdot_j.plus(Xc.get(r).times(C.get(r + 1, j)));
        // Call the DAE right-hand-side function
        List<MX> oay = (List<MX>) ffcn.invoke(Xc.get(j), Zc.get(j), Uk);
        MX Ode_j = oay.get(0), Alg_j = oay.get(1), Y_j = oay.get(2);
        // Append collocation equations
        G.add(Ode_j.times(h).minus(Xdot_j));
        // Add algebraic equations
        G.add(Alg_j);
        // Add output equations
        G.add(Y_j.minus(Yc.get(j)));
      }
      // State at the end of the interval
      MX Xk_end = Xk.times(D.get(0));
      for (int j = 0; j < deg; j++) Xk_end = Xk_end.plus(Xc.get(j).times(D.get(j + 1)));
      // New variable for the state at the end of the interval
      tk = tk + h;
      Xk = MX.sym("X" + (k + 1), X.size());
      X_all.add(Xk);
      W.add(Xk);
      Xgrid.add(tk);
      // Enforce continuity
      G.add(Xk_end.minus(Xk));
      // Output, algebraic variable at the end of the interval (only correct for Radau)
      MX Zk = Zc.get(Zc.size() - 1);
      Yk = Yc.get(Yc.size() - 1);
    }

    // Include end of control interval
    Ugrid.add(tk);

    // Add cost contribution from the final state
    MX Dx = Xk.minus(Xf_ref);
    J = J.plus(ca.dot(Dx, Dx.times(Xf_weights)));
    if (Y.size() > 0) {
      // Add contribution to the cost
      MX Dy = Yk.minus(Yf_ref);
      J = J.plus(ca.dot(Dy, Dy.times(Yf_weights)));
    }

    // Concatenate vectors
    MX Wv = ca.vertcat(W);
    MX Gv = ca.vertcat(G);

    // Allocate an NLP solver
    Map<String, MX> nlp = Map.of("x", Wv, "f", J, "g", Gv);
    Map<String, Object> opts = new LinkedHashMap<>();
    opts.put("ipopt.max_iter", 100);
    opts.put("ipopt.tol", 1e-4);
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);

    // Create mappings from (X, Z, U) -> (W) and back
    MX X_allv = ca.horzcat(X_all);
    MX Z_allv = ca.horzcat(Z_all);
    MX U_allv = ca.horzcat(U_all);
    MX Y_allv = ca.horzcat(Y_all);
    Function to_W = new Function("to_W", List.of(X_allv, Z_allv, Y_allv, U_allv), List.of(Wv),
        List.of("X", "Z", "Y", "U"), List.of("W"));
    Function from_W = new Function("from_W", List.of(Wv), List.of(X_allv, Z_allv, Y_allv, U_allv),
        List.of("W"), List.of("X", "Z", "Y", "U"));

    // Variable bounds, initial guess: Get data structures of correct size
    List<DM> o0 = (List<DM>) from_W.invoke(0);
    DM X0 = o0.get(0), Z0 = o0.get(1), Y0 = o0.get(2), U0 = o0.get(3);
    List<DM> olb = (List<DM>) from_W.invoke(-inf);
    DM lbX = olb.get(0), lbZ = olb.get(1), lbY = olb.get(2), lbU = olb.get(3);
    List<DM> oub = (List<DM>) from_W.invoke(inf);
    DM ubX = oub.get(0), ubZ = oub.get(1), ubY = oub.get(2), ubU = oub.get(3);

    // Initial guess, bounds for control trajectory
    List<Object> Uinit = dae.get(U);
    double[] Umin = dae.min(U);
    double[] Umax = dae.max(U);
    for (int i = 0; i < Uinit.size(); i++) {
      double v0 = (Double) Uinit.get(i);
      U0.set(i, new Slice(), v0);
      lbU.set(i, new Slice(), Umin[i]);
      ubU.set(i, new Slice(), Umax[i]);
    }

    // Initial guess, bounds for state trajectory
    List<Object> Xinit = dae.get(X);
    double[] Xmin = dae.min(X);
    double[] Xmax = dae.max(X);
    for (int i = 0; i < Xinit.size(); i++) {
      double v0 = (Double) Xinit.get(i);
      X0.set(i, new Slice(), v0);
      lbX.set(i, new Slice(), Xmin[i]);
      lbX.set(i, 0, v0);
      ubX.set(i, new Slice(), Xmax[i]);
      ubX.set(i, 0, v0);
    }

    // Initial guess for algebraic variable trajectory
    List<Object> Zinit = dae.get(Z);
    double[] Zmin = dae.min(Z);
    double[] Zmax = dae.max(Z);
    for (int i = 0; i < Zinit.size(); i++) {
      double v0 = (Double) Zinit.get(i);
      Z0.set(i, new Slice(), v0);
      lbZ.set(i, new Slice(), Zmin[i]);
      ubZ.set(i, new Slice(), Zmax[i]);
    }

    // Initial guess for output trajectory
    List<Object> Yinit = dae.get(Y);
    double[] Ymin = dae.min(Y);
    double[] Ymax = dae.max(Y);
    for (int i = 0; i < Yinit.size(); i++) {
      double v0 = (Double) Yinit.get(i);
      Y0.set(i, new Slice(), v0);
      lbY.set(i, new Slice(), Ymin[i]);
      ubY.set(i, new Slice(), Ymax[i]);
    }

    // Solve the NLP
    DM W0 = (DM) to_W.invoke(X0, Z0, Y0, U0);
    DM lbW = (DM) to_W.invoke(lbX, lbZ, lbY, lbU);
    DM ubW = (DM) to_W.invoke(ubX, ubZ, ubY, ubU);
    Map<String, GenericMatrix> sol = solver.invoke(Map.of("x0", W0, "lbx", lbW, "ubx", ubW, "lbg", 0, "ubg", 0));
    List<DM> opt = (List<DM>) from_W.invoke(sol.get("x"));
    double[][] Xopt = opt.get(0).full();
    double[][] Zopt = opt.get(1).full();
    double[][] Yopt = opt.get(2).full();
    double[][] Uopt = opt.get(3).full();

    // Variables to plot
    List<String> plot_vars = new ArrayList<>();
    plot_vars.addAll(X);
    plot_vars.addAll(Z);
    plot_vars.addAll(Y);
    plot_vars.addAll(U);

    // plot: one subplot per plot_vars entry: Xopt rows over Xgrid, Zopt over Zgrid, Yopt over Ygrid,
    //       Uopt (step) over Ugrid
  }
}
