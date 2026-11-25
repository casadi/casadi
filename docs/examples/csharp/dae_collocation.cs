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

public static class dae_collocation {

  static bool hasMa27() {
    SX z = SX.sym("z");
    Function s = ca.nlpsol("probe", "ipopt", new Dictionary<string, SX> { ["x"] = z, ["f"] = z * z },
        new Dictionary<string, object> { ["ipopt.linear_solver"] = "ma27", ["ipopt.print_level"] = 0L, ["print_time"] = false });
    s.call(new Dictionary<string, DM> { ["x0"] = 1 });
    return (bool)s.stats()["success"];
  }

  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }
    if (!hasMa27()) { Console.WriteLine("SKIP: ipopt linear solver ma27 (libhsl) not available"); return; }

    // -----------------------------------------------------------------------------
    // Collocation setup
    // -----------------------------------------------------------------------------
    int nicp = 1;        // Number of (intermediate) collocation points per control interval

    double xref = 0.1; // chariot reference

    double l = 1.0; //- -> crane, + -> pendulum
    double m = 1.0;
    double M = 1.0;
    double g = 9.81;
    double tf = 5.0;
    int nk = 50;
    int ndstate = 6;
    int nastate = 1;
    int ninput = 1;

    // Degree of interpolating polynomial
    int deg = 4;
    // Radau collocation points
    string cp = "radau";
    // Size of the finite elements
    double h = tf / nk / nicp;

    // Coefficients of the collocation equation
    var C = new double[deg + 1, deg + 1];
    // Coefficients of the continuity equation
    var D = new double[deg + 1];

    // Collocation point
    SX tau = SX.sym("tau");

    // All collocation time points
    double[] tau_root = ca.collocation_points(deg, cp).Prepend(0.0).ToArray();

    var T = new double[nk, deg + 1];
    for (int i = 0; i < nk; i++)
      for (int j = 0; j < deg + 1; j++)
        T[i, j] = h * (i + tau_root[j]);

    // For all collocation points: eq 10.4 or 10.17 in Biegler's book
    // Construct Lagrange polynomials to get the polynomial basis at the collocation point
    for (int j = 0; j < deg + 1; j++) {
      SX L = 1;
      for (int j2 = 0; j2 < deg + 1; j2++) {
        if (j2 != j) {
          L *= (tau - tau_root[j2]) / (tau_root[j] - tau_root[j2]);
        }
      }

      // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      var lfcn = ca.Function("lfcn", [tau], [L]);
      D[j] = (double)lfcn.Call(1.0);

      // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      var tfcn = ca.Function("tfcn", [tau], [ca.tangent(L, tau)]);
      for (int j2 = 0; j2 < deg + 1; j2++) {
        C[j, j2] = (double)tfcn.Call(tau_root[j2]);
      }
    }

    // -----------------------------------------------------------------------------
    // Model setup
    // -----------------------------------------------------------------------------
    // Declare variables (use scalar graph)
    SX t = SX.sym("t");                  // time
    SX u = SX.sym("u");                  // control
    SX xd = SX.sym("xd", ndstate);       // differential state
    SX xa = SX.sym("xa", nastate);       // algebraic state
    SX xddot = SX.sym("xdot", ndstate);  // differential state time derivative
    SX p = SX.sym("p", 0, 1);            // parameters

    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX w = SX.sym("w");

    SX dx = SX.sym("dx");
    SX dy = SX.sym("dy");
    SX dw = SX.sym("dw");


    SX res = ca.vertcat(xddot[0] - dx,
                     xddot[1] - dy,
                     xddot[2] - dw,
                     m * xddot[3] + (x - w) * xa,
                     m * xddot[4] + y * xa - g * m,
                     M * xddot[5] + (w - x) * xa + u,
                     (x - w) * (xddot[3] - xddot[5]) + y * xddot[4] + dy * dy + (dx - dw) * (dx - dw));


    xd[0] = x;
    xd[1] = y;
    xd[2] = w;
    xd[3] = dx;
    xd[4] = dy;
    xd[5] = dw;


    // System dynamics (implicit formulation)
    var ffcn = ca.Function("ffcn", [t, xddot, xd, xa, u, p], [res]);

    // Objective function
    var MayerTerm = ca.Function("mayer", [t, xd, xa, u, p],
        [(x - xref) * (x - xref) + (w - xref) * (w - xref) + dx * dx + dy * dy]);
    var LagrangeTerm = ca.Function("lagrange", [t, xd, xa, u, p],
        [(x - xref) * (x - xref) + (w - xref) * (w - xref)]);

    // Control bounds
    double[] u_min = { -2 };
    double[] u_max = { 2 };
    double[][] u_init = Enumerable.Range(0, nk * nicp * (deg + 1)).Select(_ => new double[1]).ToArray(); // needs to be specified for every time interval (even though it stays constant)

    // Differential state bounds
    // Path bounds
    double[] xD_min = { -ca.inf, -ca.inf, -ca.inf, -ca.inf, -ca.inf, -ca.inf };
    double[] xD_max = { ca.inf, ca.inf, ca.inf, ca.inf, ca.inf, ca.inf };
    // Initial bounds
    double[] xDi_min = { 0.0, l, 0.0, 0.0, 0.0, 0.0 };
    double[] xDi_max = { 0.0, l, 0.0, 0.0, 0.0, 0.0 };
    // Final bounds
    double[] xDf_min = { -ca.inf, -ca.inf, -ca.inf, -ca.inf, -ca.inf, -ca.inf };
    double[] xDf_max = { ca.inf, ca.inf, ca.inf, ca.inf, ca.inf, ca.inf };

    // Initial guess for differential states
    double[][] xD_init = Enumerable.Range(0, nk * nicp * (deg + 1))
        .Select(_ => new double[] { 0.0, l, 0.0, 0.0, 0.0, 0.0 }).ToArray(); // needs to be specified for every time interval

    // Algebraic state bounds and initial guess
    double[] xA_min = { -ca.inf };
    double[] xA_max = { ca.inf };
    double[] xAi_min = { -ca.inf };
    double[] xAi_max = { ca.inf };
    double[] xAf_min = { -ca.inf };
    double[] xAf_max = { ca.inf };
    double[][] xA_init = Enumerable.Range(0, nk * nicp * (deg + 1))
        .Select(_ => new double[] { Math.Sign(l) * 9.81 }).ToArray();

    // Parameter bounds and initial guess
    double[] p_min = { };
    double[] p_max = { };
    double[] p_init = { };




    // -----------------------------------------------------------------------------
    // Constraints setup
    // -----------------------------------------------------------------------------
    // Initial constraint
    double[] ic_min = { };
    double[] ic_max = { };
    SX ic = ca.SX();
    // ic.append();       ic_min = append(ic_min, 0.);         ic_max = append(ic_max, 0.)
    var icfcn = ca.Function("icfcn", [t, xd, xa, u, p], [ic]);
    // Path constraint
    double[] pc_min = { };
    double[] pc_max = { };
    SX pc = ca.SX();
    // pc.append();       pc_min = append(pc_min, 0.);         pc_max = append(pc_max, 0.)
    var pcfcn = ca.Function("pcfcn", [t, xd, xa, u, p], [pc]);
    // Final constraint
    double[] fc_min = { };
    double[] fc_max = { };
    SX fc = ca.SX();
    // fc.append();       fc_min = append(fc_min, 0.);         fc_max = append(fc_max, 0.)
    var fcfcn = ca.Function("fcfcn", [t, xd, xa, u, p], [fc]);

    // -----------------------------------------------------------------------------
    // NLP setup
    // -----------------------------------------------------------------------------
    // Dimensions of the problem
    int nx = (int)(xd.nnz() + xa.nnz()); // total number of states        //MODIF
    int ndiff = (int)xd.nnz();           // number of differential states //MODIF
    int nalg = (int)xa.nnz();            // number of algebraic states
    int nu = (int)u.nnz();               // number of controls
    int NP = (int)p.nnz();               // number of parameters

    // Total number of variables
    int NXD = nicp * nk * (deg + 1) * ndiff; // Collocated differential states
    int NXA = nicp * nk * deg * nalg;        // Collocated algebraic states
    int NU = nk * nu;                        // Parametrized controls
    int NXF = ndiff;                         // Final state (only the differential states)
    int NV = NXD + NXA + NU + NXF + NP;

    // NLP variable vector
    MX V = MX.sym("V", NV);

    // All variables with bounds and initial guess
    var vars_lb = new double[NV];
    var vars_ub = new double[NV];
    var vars_init = new double[NV];
    int offset = 0;

    // Get the parameters
    MX P = V[offset..(offset + NP)];
    Array.Copy(p_init, 0, vars_init, offset, NP);
    Array.Copy(p_min, 0, vars_lb, offset, NP);
    Array.Copy(p_max, 0, vars_ub, offset, NP);
    offset += NP;

    // Get collocated states and parametrized control
    var XD = new MX[nk + 1, nicp, deg + 1]; // NB: same name as above
    var XA = new MX[nk, nicp, deg]; // NB: same name as above
    var U = new MX[nk];
    int index = 0;
    for (int k = 0; k < nk; k++) {
      // Collocated states
      for (int i = 0; i < nicp; i++) {
        //
        for (int j = 0; j < deg + 1; j++) {

          // Get the expression for the state vector
          XD[k, i, j] = V[offset..(offset + ndiff)];
          if (j != 0) {
            XA[k, i, j - 1] = V[(offset + ndiff)..(offset + ndiff + nalg)];
          }
          // Add the initial condition
          index = (deg + 1) * (nicp * k + i) + j;
          if (k == 0 && j == 0 && i == 0) {
            Array.Copy(xD_init[index], 0, vars_init, offset, ndiff);

            Array.Copy(xDi_min, 0, vars_lb, offset, ndiff);
            Array.Copy(xDi_max, 0, vars_ub, offset, ndiff);
            offset += ndiff;
          } else {
            if (j != 0) {
              Array.Copy(xD_init[index].Concat(xA_init[index]).ToArray(), 0, vars_init, offset, nx);

              Array.Copy(xD_min.Concat(xA_min).ToArray(), 0, vars_lb, offset, nx);
              Array.Copy(xD_max.Concat(xA_max).ToArray(), 0, vars_ub, offset, nx);
              offset += nx;
            } else {
              Array.Copy(xD_init[index], 0, vars_init, offset, ndiff);

              Array.Copy(xD_min, 0, vars_lb, offset, ndiff);
              Array.Copy(xD_max, 0, vars_ub, offset, ndiff);
              offset += ndiff;
            }
          }
        }
      }

      // Parametrized controls
      U[k] = V[offset..(offset + nu)];
      Array.Copy(u_min, 0, vars_lb, offset, nu);
      Array.Copy(u_max, 0, vars_ub, offset, nu);
      Array.Copy(u_init[index], 0, vars_init, offset, nu);
      offset += nu;
    }

    // State at end time
    XD[nk, 0, 0] = V[offset..(offset + ndiff)];
    Array.Copy(xDf_min, 0, vars_lb, offset, ndiff);
    Array.Copy(xDf_max, 0, vars_ub, offset, ndiff);
    Array.Copy(xD_init[^1], 0, vars_init, offset, ndiff);
    offset += ndiff;
    if (offset != NV) throw new Exception("assertion failed");

    // Constraint function for the NLP
    var gl = new List<MX>();
    var lbg = new List<double>();
    var ubg = new List<double>();

    // Initial constraints
    MX ick = icfcn.Call(0.0, XD[0, 0, 0], XA[0, 0, 0], U[0], P);
    gl.Add(ick);
    lbg.AddRange(ic_min);
    ubg.AddRange(ic_max);

    // For all finite elements
    {
      for (int k = 0; k < nk; k++) {
        for (int i = 0; i < nicp; i++) {
          // For all collocation points
          for (int j = 1; j < deg + 1; j++) {
            // Get an expression for the state derivative at the collocation point
            MX xp_jk = 0;
            for (int j2 = 0; j2 < deg + 1; j2++) {
              xp_jk += C[j2, j] * XD[k, i, j2]; // get the time derivative of the differential states (eq 10.19b)
            }

            // Add collocation equations to the NLP
            MX fk = ffcn.Call(0.0, xp_jk / h, XD[k, i, j], XA[k, i, j - 1], U[k], P);
            gl.Add(fk[..ndiff]);           // impose system dynamics (for the differential states (eq 10.19b))
            lbg.AddRange(new double[ndiff]); // equality constraints
            ubg.AddRange(new double[ndiff]); // equality constraints
            gl.Add(fk[ndiff..]);           // impose system dynamics (for the algebraic states (eq 10.19b))
            lbg.AddRange(new double[nalg]); // equality constraints
            ubg.AddRange(new double[nalg]); // equality constraints

            //  Evaluate the path constraint function
            MX pck = pcfcn.Call(0.0, XD[k, i, j], XA[k, i, j - 1], U[k], P);

            gl.Add(pck);
            lbg.AddRange(pc_min);
            ubg.AddRange(pc_max);
          }

          // Get an expression for the state at the end of the finite element
          MX xf_k = 0;
          for (int j = 0; j < deg + 1; j++) {
            xf_k += D[j] * XD[k, i, j];
          }

          // Add continuity equation to NLP
          if (i == nicp - 1) {
            gl.Add(XD[k + 1, 0, 0] - xf_k);
          } else {
            gl.Add(XD[k, i + 1, 0] - xf_k);
          }

          lbg.AddRange(new double[ndiff]);
          ubg.AddRange(new double[ndiff]);
        }
      }
    }
    // Python leaves the loop variables at their last values
    {
      int k = nk - 1, i = nicp - 1, j = deg;

      // Periodicity constraints
      //   none

      // Final constraints (Const, dConst, ConstQ)
      MX fck = fcfcn.Call(0.0, XD[k, i, j], XA[k, i, j - 1], U[k], P);
      gl.Add(fck);
      lbg.AddRange(fc_min);
      ubg.AddRange(fc_max);
    }

    // Objective function of the NLP
    // Implement Mayer term
    MX Obj = 0;
    {
      int k = nk - 1, i = nicp - 1, j = deg;
      MX obj = MayerTerm.Call(0.0, XD[k, i, j], XA[k, i, j - 1], U[k], P);
      Obj += obj;
    }

    // Implement Lagrange term
    DM lDotAtTauRoot = ((DM)C).T();
    DM lAtOne = D;

    DM ldInv = ca.inv(lDotAtTauRoot[1..(deg + 1), 1..(deg + 1)]);
    DM ld0 = lDotAtTauRoot[1..(deg + 1), 0];
    MX lagrangeTerm = 0;
    for (int k = 0; k < nk; k++) {
      for (int i = 0; i < nicp; i++) {
        var lt = new List<MX>();
        for (int j = 1; j < deg + 1; j++) lt.Add(LagrangeTerm.Call(0.0, XD[k, i, j], XA[k, i, j - 1], U[k], P));
        MX dQs = h * ca.veccat(lt);
        MX Qs = ca.mtimes(ldInv, dQs);
        MX m_ = ca.mtimes(Qs.T(), lAtOne[1..(deg + 1)]);
        lagrangeTerm += m_;
      }
    }

    Obj += lagrangeTerm;

    // NLP
    var nlp = new Dictionary<string, MX> { ["x"] = V, ["f"] = Obj, ["g"] = ca.vertcat(gl) };

    // ----
    // SOLVE THE NLP
    // ----

    // NLP solver options
    var opts = new Dictionary<string, object> {
      ["expand"] = true,
      ["ipopt.max_iter"] = 1000L,
      ["ipopt.tol"] = 1e-4,
      ["ipopt.linear_solver"] = "ma27" };

    // Allocate an NLP solver
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    var arg = new Dictionary<string, DM>();

    // Initial condition
    arg["x0"] = vars_init;

    // Bounds on x
    arg["lbx"] = vars_lb;
    arg["ubx"] = vars_ub;

    // Bounds on g
    arg["lbg"] = lbg.ToArray();
    arg["ubg"] = ubg.ToArray();

    // Solve the problem
    var resl = solver.call(arg);

    // Print the optimal cost
    Console.WriteLine($"optimal cost:  {(double)resl["f"]}");

    // Retrieve the solution
    double[] v_opt = resl["x"].nonzeros().ToArray();


    // ----
    // RETRIEVE THE SOLUTION
    // ----
    var xD_opt = new double[ndiff, (deg + 1) * nicp * nk + 1];
    var xA_opt = new double[nalg, deg * nicp * nk];
    var u_opt = new double[nu, (deg + 1) * nicp * nk + 1];
    offset = 0;
    int offset2 = 0;
    int offset3 = 0;
    int offset4 = 0;

    for (int k = 0; k < nk; k++) {
      for (int i = 0; i < nicp; i++) {
        for (int j = 0; j < deg + 1; j++) {
          for (int r = 0; r < ndiff; r++) xD_opt[r, offset2] = v_opt[offset + r];
          offset2 += 1;
          offset += ndiff;
          if (j != 0) {
            for (int r = 0; r < nalg; r++) xA_opt[r, offset4] = v_opt[offset + r];
            offset4 += 1;
            offset += nalg;
          }
        }
      }
      double[] utemp = v_opt[offset..(offset + nu)];
      for (int i = 0; i < nicp; i++) {
        for (int j = 0; j < deg + 1; j++) {
          for (int r = 0; r < nu; r++) u_opt[r, offset3] = utemp[r];
          offset3 += 1;
        }
      }
      //    u_opt += v_opt[offset:offset+nu]
      offset += nu;
    }

    for (int r = 0; r < ndiff; r++) xD_opt[r, xD_opt.GetLength(1) - 1] = v_opt[offset + r];



    // The algebraic states are not defined at the first collocation point of the finite elements:
    // with the polynomials we compute them at that point
    var Da = new double[deg];
    for (int j = 1; j < deg + 1; j++) {
      // Lagrange polynomials for the algebraic states: exclude the first point
      SX La = 1;
      for (int j2 = 1; j2 < deg + 1; j2++) {
        if (j2 != j) {
          La *= (tau - tau_root[j2]) / (tau_root[j] - tau_root[j2]);
        }
      }
      var lafcn = ca.Function("lafcn", [tau], [La]);
      Da[j - 1] = (double)lafcn.Call(tau_root[0]);
    }

    int ncol = (deg + 1) * nicp * nk + 1;
    var xA_plt = new double[nalg, ncol];
    offset4 = 0;
    int offset5 = 0;
    for (int k = 0; k < nk; k++) {
      for (int i = 0; i < nicp; i++) {
        for (int j = 0; j < deg + 1; j++) {
          if (j != 0) {
            for (int r = 0; r < nalg; r++) xA_plt[r, offset5] = xA_opt[r, offset4];
            offset4 += 1;
            offset5 += 1;
          } else {
            for (int r = 0; r < nalg; r++) {
              double xa0 = 0;
              for (int jj = 0; jj < deg; jj++) {
                xa0 += Da[jj] * xA_opt[r, offset4 + jj];
              }
              xA_plt[r, offset5] = xa0;
            }
            // xA_plt[:,offset5] = xA_opt[:,offset4]
            offset5 += 1;
          }
        }
      }
    }

    for (int r = 0; r < nalg; r++) xA_plt[r, ncol - 1] = xA_plt[r, ncol - 2];




    var tgrid = new double[nk * nicp * (deg + 1) + 1];
    int n = 0;
    double last = 0;
    for (int k = 0; k < nk * nicp; k++) {
      double @base = k == 0 ? 0 : last;
      for (int j = 0; j < deg + 1; j++) tgrid[n++] = @base + tau_root[j] * h;
      last = tgrid[n - 1];
    }
    tgrid[n] = tgrid[n - 1];
    // plot: xD_opt[0], xD_opt[1], xD_opt[2] ("x", "y", "w") over tgrid
    // plot: u_opt[0] over tgrid ("Crane, inputs")
    // plot: xA_plt[0] over tgrid ("Crane, lambda")
  }
}
