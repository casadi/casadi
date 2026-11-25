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

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.IO.Compression;
using System.Linq;
using CasADi;

public static class fmu_collocation {

  // Python repr of a str -> str dict
  static string pyrepr(IEnumerable<KeyValuePair<string, object>> m) =>
    "{" + string.Join(", ", m.Select(e => $"'{e.Key}': '{e.Value}'")) + "}";

  static double[] select(List<string> names, string[] which, double value) =>
    names.Select(n => which.Contains(n) ? value : 0).ToArray();

  // Path with "." components removed, as Python's os.path.normpath
  static string normpath(string p) => string.Join("/", p.Split('/').Where(s => s != "."));

  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }

    // Precompiled FMU?
    string fmu = null; // Generate CasADi FMU
    // fmu = "casadi_symbolic_fmu"; // FMU from CasADi, use serialized CasADi expressions (default)
    // fmu = "casadi_blackbox_fmu"; // FMU from CasADi, do not use serialized CasADi expressions
    // fmu = "dymola_fmu"; // FMU from Dymola, compiled from Modelica

    DaeBuilder dae;
    if (fmu == null) {
      // Start with an empty DaeBuilder instance
      dae = ca.DaeBuilder("Crane");

      // Hard coded constants
      double mL = 1.0; // Mass of the load
      double mC = 1.0; // Mass of the crane
      double g = 9.81; // Gravity
      double length = 1.0; // Rope/pole length
      bool inverted = true; // Inverted pendulum on a cart rather than crane with hanging load

      // States
      MX xL = dae.add("xL", new Dictionary<string, object> { ["initial"] = "exact", ["start"] = 0.0,
          ["description"] = "Horizontal position of the load" });
      MX y = dae.add("y", new Dictionary<string, object> { ["initial"] = "exact", ["start"] = inverted ? length : -length,
          ["description"] = "Vertical position of the load" });
      MX xC = dae.add("xC", new Dictionary<string, object> { ["initial"] = "exact", ["start"] = 0.0,
          ["description"] = "Horizontal position of the crane" });
      MX dxL = dae.add("dxL", new Dictionary<string, object> { ["initial"] = "exact", ["start"] = 0.0,
          ["description"] = "Horizontal velocity of the load" });
      MX dy = dae.add("dy", new Dictionary<string, object> { ["initial"] = "exact", ["start"] = 0.0,
          ["description"] = "Vertical velocity of the load" });
      MX dxC = dae.add("dxC", new Dictionary<string, object> { ["initial"] = "exact", ["start"] = 0.0,
          ["description"] = "Horizontal velocity of the crane" });

      // Input
      MX u = dae.add("u", "input", new Dictionary<string, object> { ["min"] = -2.0, ["max"] = 2.0, ["start"] = 0.0,
          ["description"] = "Horizontal force applied to the crane" });

      // Algebraic variables
      MX T = dae.add("T", new Dictionary<string, object> { ["initial"] = "approx", ["start"] = inverted ? g : -g,
          ["description"] = "Tension in the rope/pole" });

      // Ordinary differential equations
      MX ddxL = (xL - xC) / length * T / mL;
      MX ddy = g - y / length * T / mL;
      MX ddxC = ((xL - xC) / length * T - u) / mC;
      dae.eq(dae.der(xL), dxL);
      dae.eq(dae.der(y), dy);
      dae.eq(dae.der(xC), dxC);
      dae.eq(dae.der(dxL), ddxL);
      dae.eq(dae.der(dy), ddy);
      dae.eq(dae.der(dxC), ddxC);

      // Original algebraic equation
      MX res = ca.pow(y, 2) + ca.pow(xL - xC, 2) - Math.Pow(length, 2);
      Console.WriteLine($"Original residual equation: {res}");

      // Index reduction: Differentiate with respect to time
      // In CasADi, the most efficient way to do index reduction it to use forward
      // mode AD, e.g. us using Jacobian-times-vector products (jtimes):
      // Using the chain rule we have: d(res)/dt = d(res)/dx * dx/dt
      MX res1 = ca.jtimes(res, ca.vertcat(dae.inputs("x")), ca.vertcat(dae.outputs("ode")));
      MX res2 = ca.jtimes(res1, ca.vertcat(dae.inputs("x")), ca.vertcat(dae.outputs("ode")));
      Console.WriteLine($"Residual equation after index reduction: {res2}");

      // Add algebraic equations
      dae.eq(0, res2);

      // Simulation duration
      dae.set_start_time(0);
      dae.set_stop_time(5);

      // Generate FMU
      var fmu_files = dae.export_fmu().OrderBy(e => e.Key, StringComparer.Ordinal).ToList();
      Console.WriteLine($"Generated files: {pyrepr(fmu_files)}");

      // Compile DLL
      string fmi_headers = "../../../external_packages/FMI-Standard-3.0/headers";
      string cfiles = string.Join(" ", fmu_files.Select(e => e.Key).Where(f => f.EndsWith(".c")));
      // Compile sources (in Linux)
      if (Path.DirectorySeparatorChar == '/') {
        try {
          string sofile = dae.name() + ".so";
          Process.Start("sh", ["-c", $"gcc --shared -fPIC -I{fmi_headers} {cfiles} -o {sofile}"]).WaitForExit();
          Console.WriteLine($"Compiled {sofile}");
          fmu_files.Add(new KeyValuePair<string, object>(sofile, "binaries/x86_64-linux/" + sofile));
        } catch (Exception e) {
          Console.WriteLine($"Failed to create FMU binary: {e.Message}");
        }
      }

      // Package into an FMU
      string fmuname = dae.name() + ".fmu";
      File.Delete(fmuname);  // overwrite, as Python's zipfile.ZipFile(fmuname, "w")
      using (ZipArchive fmufile = ZipFile.Open(fmuname, ZipArchiveMode.Create)) {
        foreach (var e in fmu_files) {
          fmufile.CreateEntryFromFile(e.Key, normpath((string)e.Value));
          File.Delete(e.Key);
        }
      }
      Console.WriteLine($"Created FMU: {fmuname}");
    } else if (fmu == "casadi_symbolic_fmu") {
      // Use FMU generated from CasADi using the script above: Standard FMU + CasADi expressions + LS-DAE
      dae = ca.DaeBuilder("Crane", "./Crane.fmu");
    } else if (fmu == "casadi_blackbox_fmu") {
      // Use FMU generated from CasADi using the script above: Standard FMU + LS-DAE
      dae = ca.DaeBuilder("Crane", "./Crane.fmu", new Dictionary<string, object> { ["enable_ls_serialization"] = false });
    } else if (fmu == "dymola_fmu") {
      // Use FMU compiled from Modelica using Dymola
      dae = ca.DaeBuilder("crane", "../../../test/data/Crane.fmu");
    } else {
      throw new ArgumentException($"Unknown FMU option: {fmu}");
    }

    // Display the model
    dae.disp(true);

    // DAE right-hand-side function
    Function ffcn = dae.create("ffcn", ["x", "z", "u"], ["ode", "alg", "y"]);

    // Variable names
    List<string> U = dae.u();
    List<string> X = dae.x();
    List<string> Z = dae.z();
    List<string> Y = dae.y();

    // Refence location for crane and load
    double xref = 0.1;

    // Quadratic cost at the beginning of control intervals
    DM X_ref = select(X, ["xL", "xC"], xref);
    DM Y_ref = select(Y, ["xL", "xC"], xref);
    DM X_weights = select(X, ["xL", "xC"], 1.0);
    DM Y_weights = select(Y, ["xL", "xC"], 1.0);

    // Quadratic cost at final time
    DM Xf_ref = select(X, ["xL", "xC"], xref);
    DM Yf_ref = select(Y, ["xL", "xC"], xref);
    DM Xf_weights = select(X, ["xL", "xC", "dxL", "dy"], 1.0);
    DM Yf_weights = select(Y, ["xL", "xC", "dxL", "dy"], 1.0);

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
    List<double> tau_root = ca.collocation_points(deg, "radau");

    // We can query the coefficient for the interpolating polynomials using a helper function
    // C: Coefficients for calculating state derivatives
    // D: Coefficients for calculating state at the end of the interval
    // B: Coefficients for calculating a quadrature at the end of the interval
    var (C, D, B) = ca.collocation_coeff(tau_root);

    // We will construct an NLP incrementally, starting with no decision variables,
    // and constraints, and zero objective
    var W = new List<MX>();
    var G = new List<MX>();
    MX J = 0;

    // Solution trajectories
    var U_all = new List<MX>();
    var X_all = new List<MX>();
    var Z_all = new List<MX>();
    var Y_all = new List<MX>();

    // Time grids, for solution
    var Ugrid = new List<double>();
    var Xgrid = new List<double>();
    var Zgrid = new List<double>();
    var Ygrid = new List<double>();

    // Current time
    double tk = 0;

    // Add a variable for the initial state
    MX Xk = MX.sym("X0", X.Count);
    X_all.Add(Xk);
    W.Add(Xk);
    Xgrid.Add(tk);
    MX Yk = null;
    // Loop over all control intervals, constructing the NLP and x, u trajectories
    for (int k = 0; k < nk; k++) {
      // Add contribution to the cost
      MX Dx = Xk - X_ref;
      J += ca.dot(Dx, Dx * X_weights);
      if (k > 0 && Y.Count > 0) {
        // Add contribution to the cost
        MX Dy = Yk - Y_ref;
        J += ca.dot(Dy, Dy * Y_weights);
      }
      // Symbolic expression for the control for the interval
      MX Uk = MX.sym($"u{k}", U.Count);
      U_all.Add(Uk);
      W.Add(Uk);
      Ugrid.Add(tk);
      // Declare variables at each collocation point
      var Xc = new List<MX>();
      var Zc = new List<MX>();
      var Yc = new List<MX>();
      for (int j = 1; j < deg + 1; j++) {
        // Differential state at collocation points
        MX Xkj = MX.sym($"X{k}_{j}", X.Count);
        Xc.Add(Xkj);
        X_all.Add(Xkj);
        W.Add(Xkj);
        Xgrid.Add(tk + tau_root[j - 1] * h);
        // Algebraic variables at collocation points
        MX Zkj = MX.sym($"Z{k}_{j}", Z.Count);
        Zc.Add(Zkj);
        Z_all.Add(Zkj);
        W.Add(Zkj);
        Zgrid.Add(tk + tau_root[j - 1] * h);
        // Output variables at collocation points
        MX Ykj = MX.sym($"Y{k}_{j}", Y.Count);
        Yc.Add(Ykj);
        Y_all.Add(Ykj);
        W.Add(Ykj);
        Ygrid.Add(tk + tau_root[j - 1] * h);
      }
      // Collect the collocation equations
      for (int j = 0; j < deg; j++) {
        // Expression for the state derivative at the collocation point
        MX Xdot_j = C[0, j] * Xk;
        for (int r = 0; r < deg; r++) Xdot_j += C[r + 1, j] * Xc[r];
        // Call the DAE right-hand-side function
        List<MX> oay = ffcn.call([Xc[j], Zc[j], Uk]);
        MX Ode_j = oay[0], Alg_j = oay[1], Y_j = oay[2];
        // Append collocation equations
        G.Add(h * Ode_j - Xdot_j);
        // Add algebraic equations
        G.Add(Alg_j);
        // Add output equations
        G.Add(Y_j - Yc[j]);
      }
      // State at the end of the interval
      MX Xk_end = D[0] * Xk;
      for (int j = 0; j < deg; j++) Xk_end += D[j + 1] * Xc[j];
      // New variable for the state at the end of the interval
      tk = tk + h;
      Xk = MX.sym($"X{k + 1}", X.Count);
      X_all.Add(Xk);
      W.Add(Xk);
      Xgrid.Add(tk);
      // Enforce continuity
      G.Add(Xk_end - Xk);
      // Output, algebraic variable at the end of the interval (only correct for Radau)
      MX Zk = Zc[^1];
      Yk = Yc[^1];
    }

    // Include end of control interval
    Ugrid.Add(tk);

    // Add cost contribution from the final state
    {
      MX Dx = Xk - Xf_ref;
      J += ca.dot(Dx, Dx * Xf_weights);
      if (Y.Count > 0) {
        // Add contribution to the cost
        MX Dy = Yk - Yf_ref;
        J += ca.dot(Dy, Dy * Yf_weights);
      }
    }

    // Concatenate vectors
    MX Wv = ca.vertcat(W);
    MX Gv = ca.vertcat(G);

    // Allocate an NLP solver
    var nlp = new Dictionary<string, MX> { ["x"] = Wv, ["f"] = J, ["g"] = Gv };
    var opts = new Dictionary<string, object> { ["ipopt.max_iter"] = 100L, ["ipopt.tol"] = 1e-4 };
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);

    // Create mappings from (X, Z, U) -> (W) and back
    MX X_allv = ca.horzcat(X_all);
    MX Z_allv = ca.horzcat(Z_all);
    MX U_allv = ca.horzcat(U_all);
    MX Y_allv = ca.horzcat(Y_all);
    var to_W = ca.Function("to_W", [X_allv, Z_allv, Y_allv, U_allv], [Wv],
        ["X", "Z", "Y", "U"], ["W"]);
    var from_W = ca.Function("from_W", [Wv], [X_allv, Z_allv, Y_allv, U_allv],
        ["W"], ["X", "Z", "Y", "U"]);

    // Variable bounds, initial guess: Get data structures of correct size
    List<DM> o0 = from_W.call(new DM[] { 0 });
    DM X0 = o0[0], Z0 = o0[1], Y0 = o0[2], U0 = o0[3];
    List<DM> olb = from_W.call(new DM[] { -ca.inf });
    DM lbX = olb[0], lbZ = olb[1], lbY = olb[2], lbU = olb[3];
    List<DM> oub = from_W.call(new DM[] { ca.inf });
    DM ubX = oub[0], ubZ = oub[1], ubY = oub[2], ubU = oub[3];

    // Initial guess, bounds for control trajectory
    List<object> Uinit = dae.get(U);
    List<double> Umin = dae.min(U);
    List<double> Umax = dae.max(U);
    for (int i = 0; i < Uinit.Count; i++) {
      double v0 = (double)Uinit[i];
      U0[i, ..] = v0;
      lbU[i, ..] = Umin[i];
      ubU[i, ..] = Umax[i];
    }

    // Initial guess, bounds for state trajectory
    List<object> Xinit = dae.get(X);
    List<double> Xmin = dae.min(X);
    List<double> Xmax = dae.max(X);
    for (int i = 0; i < Xinit.Count; i++) {
      double v0 = (double)Xinit[i];
      X0[i, ..] = v0;
      lbX[i, ..] = Xmin[i];
      lbX[i, 0] = v0;
      ubX[i, ..] = Xmax[i];
      ubX[i, 0] = v0;
    }

    // Initial guess for algebraic variable trajectory
    List<object> Zinit = dae.get(Z);
    List<double> Zmin = dae.min(Z);
    List<double> Zmax = dae.max(Z);
    for (int i = 0; i < Zinit.Count; i++) {
      double v0 = (double)Zinit[i];
      Z0[i, ..] = v0;
      lbZ[i, ..] = Zmin[i];
      ubZ[i, ..] = Zmax[i];
    }

    // Initial guess for output trajectory
    List<object> Yinit = dae.get(Y);
    List<double> Ymin = dae.min(Y);
    List<double> Ymax = dae.max(Y);
    for (int i = 0; i < Yinit.Count; i++) {
      double v0 = (double)Yinit[i];
      Y0[i, ..] = v0;
      lbY[i, ..] = Ymin[i];
      ubY[i, ..] = Ymax[i];
    }

    // Solve the NLP
    DM W0 = to_W.Call(X0, Z0, Y0, U0);
    DM lbW = to_W.Call(lbX, lbZ, lbY, lbU);
    DM ubW = to_W.Call(ubX, ubZ, ubY, ubU);
    var sol = solver.call(new Dictionary<string, DM> { ["x0"] = W0, ["lbx"] = lbW, ["ubx"] = ubW, ["lbg"] = 0, ["ubg"] = 0 });
    List<DM> opt = from_W.call(new[] { sol["x"] });
    double[,] Xopt = opt[0].ToArray();
    double[,] Zopt = opt[1].ToArray();
    double[,] Yopt = opt[2].ToArray();
    double[,] Uopt = opt[3].ToArray();

    // Variables to plot
    var plot_vars = new List<string>();
    plot_vars.AddRange(X);
    plot_vars.AddRange(Z);
    plot_vars.AddRange(Y);
    plot_vars.AddRange(U);

    // plot: one subplot per plot_vars entry: Xopt rows over Xgrid, Zopt over Zgrid, Yopt over Ygrid,
    //       Uopt (step) over Ugrid
  }
}
