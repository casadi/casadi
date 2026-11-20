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

// @author: Mario Zanon and Sebastien Gross, KU Leuven 2012

import org.casadi.*;
import java.util.*;

public class dae_collocation {

  static boolean hasMa27() {
    SX z = SX.sym("z");
    Function s = ca.nlpsol("probe", "ipopt", Map.of("x", z, "f", z.times(z)),
        Map.of("ipopt.linear_solver", "ma27", "ipopt.print_level", 0, "print_time", false));
    s.invoke(Map.of("x0", 1));
    return (Boolean) s.stats().get("success");
  }

  static double[] append(double[] a, double[] b) {
    double[] r = Arrays.copyOf(a, a.length + b.length);
    System.arraycopy(b, 0, r, a.length, b.length);
    return r;
  }

  static double[] concat(List<double[]> v) {
    return v.stream().flatMapToDouble(Arrays::stream).toArray();
  }

  static double[] filled(int n, double v) {
    double[] r = new double[n];
    Arrays.fill(r, v);
    return r;
  }

  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }
    if (!hasMa27()) { System.out.println("SKIP: ipopt linear solver ma27 (libhsl) not available"); return; }
    double inf = Double.POSITIVE_INFINITY;

    // -----------------------------------------------------------------------------
    // Collocation setup
    // -----------------------------------------------------------------------------
    int nicp = 1;        // Number of (intermediate) collocation points per control interval

    double xref = 0.1; // chariot reference

    double l = 1.; //- -> crane, + -> pendulum
    double m = 1.;
    double M = 1.;
    double g = 9.81;
    double tf = 5.0;
    int nk = 50;
    int ndstate = 6;
    int nastate = 1;
    int ninput = 1;

    // Degree of interpolating polynomial
    int deg = 4;
    // Radau collocation points
    String cp = "radau";
    // Size of the finite elements
    double h = tf / nk / nicp;

    // Coefficients of the collocation equation
    double[][] C = new double[deg + 1][deg + 1];
    // Coefficients of the continuity equation
    double[] D = new double[deg + 1];

    // Collocation point
    SX tau = SX.sym("tau");

    // All collocation time points
    double[] tau_root = new double[deg + 1];
    System.arraycopy(ca.collocation_points(deg, cp), 0, tau_root, 1, deg);

    double[][] T = new double[nk][deg + 1];
    for (int i = 0; i < nk; i++)
      for (int j = 0; j < deg + 1; j++)
        T[i][j] = h * (i + tau_root[j]);

    // For all collocation points: eq 10.4 or 10.17 in Biegler's book
    // Construct Lagrange polynomials to get the polynomial basis at the collocation point
    for (int j = 0; j < deg + 1; j++) {
      SX L = new SX(1);
      for (int j2 = 0; j2 < deg + 1; j2++) {
        if (j2 != j) {
          L = L.times(tau.minus(tau_root[j2]).div(tau_root[j] - tau_root[j2]));
        }
      }

      // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      Function lfcn = new Function("lfcn", List.of(tau), List.of(L));
      D[j] = ((DM) lfcn.invoke(1.0)).to_double();

      // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      Function tfcn = new Function("tfcn", List.of(tau), List.of(ca.tangent(L, tau)));
      for (int j2 = 0; j2 < deg + 1; j2++) {
        C[j][j2] = ((DM) tfcn.invoke(tau_root[j2])).to_double();
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


    SX res = ca.vertcat(List.of(xddot.get(0).minus(dx),
        xddot.get(1).minus(dy),
        xddot.get(2).minus(dw),
        xddot.get(3).times(m).plus(x.minus(w).times(xa)),
        xddot.get(4).times(m).plus(y.times(xa)).minus(g * m),
        xddot.get(5).times(M).plus(w.minus(x).times(xa)).plus(u),
        x.minus(w).times(xddot.get(3).minus(xddot.get(5))).plus(y.times(xddot.get(4))).plus(dy.times(dy))
            .plus(dx.minus(dw).times(dx.minus(dw)))));


    xd.set(0, x);
    xd.set(1, y);
    xd.set(2, w);
    xd.set(3, dx);
    xd.set(4, dy);
    xd.set(5, dw);


    // System dynamics (implicit formulation)
    Function ffcn = new Function("ffcn", List.of(t, xddot, xd, xa, u, p), List.of(res));

    // Objective function
    Function MayerTerm = new Function("mayer", List.of(t, xd, xa, u, p), List.of(x.minus(xref).times(x.minus(xref))
        .plus(w.minus(xref).times(w.minus(xref))).plus(dx.times(dx)).plus(dy.times(dy))));
    Function LagrangeTerm = new Function("lagrange", List.of(t, xd, xa, u, p), List.of(x.minus(xref).times(x.minus(xref))
        .plus(w.minus(xref).times(w.minus(xref)))));

    // Control bounds
    double[] u_min = {-2};
    double[] u_max = {2};
    double[][] u_init = new double[nk * nicp * (deg + 1)][1]; // needs to be specified for every time interval (even though it stays constant)

    // Differential state bounds
    // Path bounds
    double[] xD_min = {-inf, -inf, -inf, -inf, -inf, -inf};
    double[] xD_max = { inf,  inf,  inf,  inf,  inf,  inf};
    // Initial bounds
    double[] xDi_min = {0.0, l, 0.0, 0.0, 0.0, 0.0};
    double[] xDi_max = {0.0, l, 0.0, 0.0, 0.0, 0.0};
    // Final bounds
    double[] xDf_min = {-inf, -inf, -inf, -inf, -inf, -inf};
    double[] xDf_max = { inf,  inf,  inf,  inf,  inf,  inf};

    // Initial guess for differential states
    double[][] xD_init = new double[nk * nicp * (deg + 1)][]; // needs to be specified for every time interval
    for (int r = 0; r < xD_init.length; r++) xD_init[r] = new double[] {0.0, l, 0.0, 0.0, 0.0, 0.0};

    // Algebraic state bounds and initial guess
    double[] xA_min = {-inf};
    double[] xA_max = { inf};
    double[] xAi_min = {-inf};
    double[] xAi_max = { inf};
    double[] xAf_min = {-inf};
    double[] xAf_max = { inf};
    double[][] xA_init = new double[nk * nicp * (deg + 1)][];
    for (int r = 0; r < xA_init.length; r++) xA_init[r] = new double[] {Math.signum(l) * 9.81};

    // Parameter bounds and initial guess
    double[] p_min = {};
    double[] p_max = {};
    double[] p_init = {};




    // -----------------------------------------------------------------------------
    // Constraints setup
    // -----------------------------------------------------------------------------
    // Initial constraint
    double[] ic_min = {};
    double[] ic_max = {};
    SX ic = new SX();
    // ic.append();       ic_min = append(ic_min, 0.);         ic_max = append(ic_max, 0.)
    Function icfcn = new Function("icfcn", List.of(t, xd, xa, u, p), List.of(ic));
    // Path constraint
    double[] pc_min = {};
    double[] pc_max = {};
    SX pc = new SX();
    // pc.append();       pc_min = append(pc_min, 0.);         pc_max = append(pc_max, 0.)
    Function pcfcn = new Function("pcfcn", List.of(t, xd, xa, u, p), List.of(pc));
    // Final constraint
    double[] fc_min = {};
    double[] fc_max = {};
    SX fc = new SX();
    // fc.append();       fc_min = append(fc_min, 0.);         fc_max = append(fc_max, 0.)
    Function fcfcn = new Function("fcfcn", List.of(t, xd, xa, u, p), List.of(fc));

    // -----------------------------------------------------------------------------
    // NLP setup
    // -----------------------------------------------------------------------------
    // Dimensions of the problem
    int nx = (int) (xd.nnz() + xa.nnz()); // total number of states        //MODIF
    int ndiff = (int) xd.nnz();           // number of differential states //MODIF
    int nalg = (int) xa.nnz();            // number of algebraic states
    int nu = (int) u.nnz();               // number of controls
    int NP = (int) p.nnz();               // number of parameters

    // Total number of variables
    int NXD = nicp * nk * (deg + 1) * ndiff; // Collocated differential states
    int NXA = nicp * nk * deg * nalg;        // Collocated algebraic states
    int NU = nk * nu;                        // Parametrized controls
    int NXF = ndiff;                         // Final state (only the differential states)
    int NV = NXD + NXA + NU + NXF + NP;

    // NLP variable vector
    MX V = MX.sym("V", NV);

    // All variables with bounds and initial guess
    double[] vars_lb = new double[NV];
    double[] vars_ub = new double[NV];
    double[] vars_init = new double[NV];
    int offset = 0;

    // Get the parameters
    MX P = V.get(new Slice(offset, offset + NP));
    System.arraycopy(p_init, 0, vars_init, offset, NP);
    System.arraycopy(p_min, 0, vars_lb, offset, NP);
    System.arraycopy(p_max, 0, vars_ub, offset, NP);
    offset += NP;

    // Get collocated states and parametrized control
    MX[][][] XD = new MX[nk + 1][nicp][deg + 1]; // NB: same name as above
    MX[][][] XA = new MX[nk][nicp][deg]; // NB: same name as above
    MX[] U = new MX[nk];
    int index = 0;
    for (int k = 0; k < nk; k++) {
      // Collocated states
      for (int i = 0; i < nicp; i++) {
        //
        for (int j = 0; j < deg + 1; j++) {

          // Get the expression for the state vector
          XD[k][i][j] = V.get(new Slice(offset, offset + ndiff));
          if (j != 0) {
            XA[k][i][j - 1] = V.get(new Slice(offset + ndiff, offset + ndiff + nalg));
          }
          // Add the initial condition
          index = (deg + 1) * (nicp * k + i) + j;
          if (k == 0 && j == 0 && i == 0) {
            System.arraycopy(xD_init[index], 0, vars_init, offset, ndiff);

            System.arraycopy(xDi_min, 0, vars_lb, offset, ndiff);
            System.arraycopy(xDi_max, 0, vars_ub, offset, ndiff);
            offset += ndiff;
          } else {
            if (j != 0) {
              System.arraycopy(append(xD_init[index], xA_init[index]), 0, vars_init, offset, nx);

              System.arraycopy(append(xD_min, xA_min), 0, vars_lb, offset, nx);
              System.arraycopy(append(xD_max, xA_max), 0, vars_ub, offset, nx);
              offset += nx;
            } else {
              System.arraycopy(xD_init[index], 0, vars_init, offset, ndiff);

              System.arraycopy(xD_min, 0, vars_lb, offset, ndiff);
              System.arraycopy(xD_max, 0, vars_ub, offset, ndiff);
              offset += ndiff;
            }
          }
        }
      }

      // Parametrized controls
      U[k] = V.get(new Slice(offset, offset + nu));
      System.arraycopy(u_min, 0, vars_lb, offset, nu);
      System.arraycopy(u_max, 0, vars_ub, offset, nu);
      System.arraycopy(u_init[index], 0, vars_init, offset, nu);
      offset += nu;
    }

    // State at end time
    XD[nk][0][0] = V.get(new Slice(offset, offset + ndiff));
    System.arraycopy(xDf_min, 0, vars_lb, offset, ndiff);
    System.arraycopy(xDf_max, 0, vars_ub, offset, ndiff);
    System.arraycopy(xD_init[xD_init.length - 1], 0, vars_init, offset, ndiff);
    offset += ndiff;
    if (offset != NV) throw new AssertionError();

    // Constraint function for the NLP
    List<MX> gl = new ArrayList<>();
    List<double[]> lbg = new ArrayList<>();
    List<double[]> ubg = new ArrayList<>();

    // Initial constraints
    MX ick = (MX) icfcn.invoke(0., XD[0][0][0], XA[0][0][0], U[0], P);
    gl.add(ick);
    lbg.add(ic_min);
    ubg.add(ic_max);

    // For all finite elements
    int k = 0, i = 0, j = 0;
    for (k = 0; k < nk; k++) {
      for (i = 0; i < nicp; i++) {
        // For all collocation points
        for (j = 1; j < deg + 1; j++) {
          // Get an expression for the state derivative at the collocation point
          MX xp_jk = new MX(0);
          for (int j2 = 0; j2 < deg + 1; j2++) {
            xp_jk = xp_jk.plus(XD[k][i][j2].times(C[j2][j])); // get the time derivative of the differential states (eq 10.19b)
          }

          // Add collocation equations to the NLP
          MX fk = (MX) ffcn.invoke(0., xp_jk.div(h), XD[k][i][j], XA[k][i][j - 1], U[k], P);
          gl.add(fk.get(new Slice(0, ndiff)));          // impose system dynamics (for the differential states (eq 10.19b))
          lbg.add(new double[ndiff]); // equality constraints
          ubg.add(new double[ndiff]); // equality constraints
          gl.add(fk.get(new Slice(ndiff, fk.size1()))); // impose system dynamics (for the algebraic states (eq 10.19b))
          lbg.add(new double[nalg]); // equality constraints
          ubg.add(new double[nalg]); // equality constraints

          //  Evaluate the path constraint function
          MX pck = (MX) pcfcn.invoke(0., XD[k][i][j], XA[k][i][j - 1], U[k], P);

          gl.add(pck);
          lbg.add(pc_min);
          ubg.add(pc_max);
        }

        // Get an expression for the state at the end of the finite element
        MX xf_k = new MX(0);
        for (j = 0; j < deg + 1; j++) {
          xf_k = xf_k.plus(XD[k][i][j].times(D[j]));
        }

        // Add continuity equation to NLP
        if (i == nicp - 1) {
          gl.add(XD[k + 1][0][0].minus(xf_k));
        } else {
          gl.add(XD[k][i + 1][0].minus(xf_k));
        }

        lbg.add(new double[ndiff]);
        ubg.add(new double[ndiff]);
      }
    }
    // Python leaves the loop variables at their last values
    k = nk - 1; i = nicp - 1; j = deg;

    // Periodicity constraints
    //   none

    // Final constraints (Const, dConst, ConstQ)
    MX fck = (MX) fcfcn.invoke(0., XD[k][i][j], XA[k][i][j - 1], U[k], P);
    gl.add(fck);
    lbg.add(fc_min);
    ubg.add(fc_max);

    // Objective function of the NLP
    // Implement Mayer term
    MX Obj = new MX(0);
    MX obj = (MX) MayerTerm.invoke(0., XD[k][i][j], XA[k][i][j - 1], U[k], P);
    Obj = Obj.plus(obj);

    // Implement Lagrange term
    DM lDotAtTauRoot = new DM(C).T();
    DM lAtOne = new DM(D);

    DM ldInv = ca.inv(lDotAtTauRoot.get(new Slice(1, deg + 1), new Slice(1, deg + 1)));
    DM ld0 = lDotAtTauRoot.get(new Slice(1, deg + 1), 0);
    MX lagrangeTerm = new MX(0);
    for (k = 0; k < nk; k++) {
      for (i = 0; i < nicp; i++) {
        List<MX> lt = new ArrayList<>();
        for (j = 1; j < deg + 1; j++) lt.add((MX) LagrangeTerm.invoke(0., XD[k][i][j], XA[k][i][j - 1], U[k], P));
        MX dQs = ca.veccat(lt).times(h);
        MX Qs = ca.mtimes(ldInv, dQs);
        MX m_ = ca.mtimes(Qs.T(), lAtOne.get(new Slice(1, deg + 1)));
        lagrangeTerm = lagrangeTerm.plus(m_);
      }
    }

    Obj = Obj.plus(lagrangeTerm);

    // NLP
    Map<String, MX> nlp = Map.of("x", V, "f", Obj, "g", ca.vertcat(gl));

    // ----
    // SOLVE THE NLP
    // ----

    // NLP solver options
    Map<String, Object> opts = new LinkedHashMap<>();
    opts.put("expand", true);
    opts.put("ipopt.max_iter", 1000);
    opts.put("ipopt.tol", 1e-4);
    opts.put("ipopt.linear_solver", "ma27");

    // Allocate an NLP solver
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    Map<String, Object> arg = new LinkedHashMap<>();

    // Initial condition
    arg.put("x0", vars_init);

    // Bounds on x
    arg.put("lbx", vars_lb);
    arg.put("ubx", vars_ub);

    // Bounds on g
    arg.put("lbg", concat(lbg));
    arg.put("ubg", concat(ubg));

    // Solve the problem
    Map<String, GenericMatrix> resl = solver.invoke(arg);

    // Print the optimal cost
    System.out.println("optimal cost:  " + ((DM) resl.get("f")).to_double());

    // Retrieve the solution
    double[] v_opt = ((DM) resl.get("x")).nonzeros();


    // ----
    // RETRIEVE THE SOLUTION
    // ----
    double[][] xD_opt = new double[ndiff][(deg + 1) * nicp * nk + 1];
    double[][] xA_opt = new double[nalg][deg * nicp * nk];
    double[][] u_opt = new double[nu][(deg + 1) * nicp * nk + 1];
    offset = 0;
    int offset2 = 0;
    int offset3 = 0;
    int offset4 = 0;

    for (k = 0; k < nk; k++) {
      for (i = 0; i < nicp; i++) {
        for (j = 0; j < deg + 1; j++) {
          for (int r = 0; r < ndiff; r++) xD_opt[r][offset2] = v_opt[offset + r];
          offset2 += 1;
          offset += ndiff;
          if (j != 0) {
            for (int r = 0; r < nalg; r++) xA_opt[r][offset4] = v_opt[offset + r];
            offset4 += 1;
            offset += nalg;
          }
        }
      }
      double[] utemp = Arrays.copyOfRange(v_opt, offset, offset + nu);
      for (i = 0; i < nicp; i++) {
        for (j = 0; j < deg + 1; j++) {
          for (int r = 0; r < nu; r++) u_opt[r][offset3] = utemp[r];
          offset3 += 1;
        }
      }
      //    u_opt += v_opt[offset:offset+nu]
      offset += nu;
    }

    for (int r = 0; r < ndiff; r++) xD_opt[r][xD_opt[r].length - 1] = v_opt[offset + r];



    // The algebraic states are not defined at the first collocation point of the finite elements:
    // with the polynomials we compute them at that point
    double[] Da = new double[deg];
    for (j = 1; j < deg + 1; j++) {
      // Lagrange polynomials for the algebraic states: exclude the first point
      SX La = new SX(1);
      for (int j2 = 1; j2 < deg + 1; j2++) {
        if (j2 != j) {
          La = La.times(tau.minus(tau_root[j2]).div(tau_root[j] - tau_root[j2]));
        }
      }
      Function lafcn = new Function("lafcn", List.of(tau), List.of(La));
      Da[j - 1] = ((DM) lafcn.invoke(tau_root[0])).to_double();
    }

    int ncol = (deg + 1) * nicp * nk + 1;
    double[][] xA_plt = new double[nalg][ncol];
    offset4 = 0;
    int offset5 = 0;
    for (k = 0; k < nk; k++) {
      for (i = 0; i < nicp; i++) {
        for (j = 0; j < deg + 1; j++) {
          if (j != 0) {
            for (int r = 0; r < nalg; r++) xA_plt[r][offset5] = xA_opt[r][offset4];
            offset4 += 1;
            offset5 += 1;
          } else {
            for (int r = 0; r < nalg; r++) {
              double xa0 = 0;
              for (int jj = 0; jj < deg; jj++) {
                xa0 += Da[jj] * xA_opt[r][offset4 + jj];
              }
              xA_plt[r][offset5] = xa0;
            }
            // xA_plt[:,offset5] = xA_opt[:,offset4]
            offset5 += 1;
          }
        }
      }
    }

    for (int r = 0; r < nalg; r++) xA_plt[r][ncol - 1] = xA_plt[r][ncol - 2];




    double[] tgrid = new double[nk * nicp * (deg + 1) + 1];
    int n = 0;
    double last = 0;
    for (k = 0; k < nk * nicp; k++) {
      double base = k == 0 ? 0 : last;
      for (j = 0; j < deg + 1; j++) tgrid[n++] = base + tau_root[j] * h;
      last = tgrid[n - 1];
    }
    tgrid[n] = tgrid[n - 1];
    // plot: xD_opt[0], xD_opt[1], xD_opt[2] ("x", "y", "w") over tgrid
    // plot: u_opt[0] over tgrid ("Crane, inputs")
    // plot: xA_plt[0] over tgrid ("Crane, lambda")
  }
}
