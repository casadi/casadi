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
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.SX;
import org.casadi.ca;

public class rocket_ipopt {
  public static void main(String[] args) throws FileNotFoundException {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }
    System.out.println("program started");

    // Dimensions
    int nu = 20;  // Number of control segments
    int nj = 100; // Number of integration steps per control segment

    // optimization variable
    SX u = SX.sym("u", nu); // control

    SX s_0 = new SX(0); // initial position
    SX v_0 = new SX(0); // initial speed
    SX m_0 = new SX(1); // initial mass

    SX dt = new SX(10.0 / (nj * nu)); // time step
    SX alpha = new SX(0.05); // friction
    SX beta = new SX(0.1); // fuel consumption rate

    // Trajectory
    List<SX> s_k = new ArrayList<>(), v_k = new ArrayList<>(), m_k = new ArrayList<>();

    // Integrate over the interval with Euler forward
    SX s = s_0, v = v_0, m = m_0;
    for (int k = 0; k < nu; ++k) {
      for (int j = 0; j < nj; ++j) {
        s = s.plus(dt.times(v));
        v = v.plus(dt.div(m).times(u.get(k).minus(alpha.times(v).times(v))));
        m = m.plus(dt.unaryMinus().times(beta).times(u.get(k)).times(u.get(k)));
      }
      s_k.add(s);
      v_k.add(v);
      m_k.add(m);
    }
    SX s_all = ca.vertcat(s_k), v_all = ca.vertcat(v_k), m_all = ca.vertcat(m_k);

    // Objective function
    SX f = ca.dot(u, u);

    // Terminal constraints
    SX g = ca.vertcat(s, v, v_all);

    // Create the NLP
    Map<String, SX> nlp = Map.of("x", u, "f", f, "g", g);

    // Allocate an NLP solver and buffers
    Function solver = ca.nlpsol("solver", "ipopt", nlp);

    // Bounds on g
    double[] gmin = Arrays.copyOf(new double[] {10, 0}, 2 + nu);
    double[] gmax = Arrays.copyOf(new double[] {10, 0}, 2 + nu);
    Arrays.fill(gmin, 2, 2 + nu, Double.NEGATIVE_INFINITY);
    Arrays.fill(gmax, 2, 2 + nu, 1.1);

    // Solve the problem
    Map<String, DM> arg = Map.of("lbx", new DM(-10),
                                 "ubx", new DM(10),
                                 "x0", new DM(0.4),
                                 "lbg", new DM(gmin),
                                 "ubg", new DM(gmax));
    Map<String, DM> res = solver.call(arg);

    // Print the optimal cost
    double cost = res.get("f").to_double();
    System.out.println("optimal cost: " + str(cost));

    // Print the optimal solution
    double[] uopt = res.get("x").nonzeros();
    System.out.println("optimal control: " + str(uopt));

    // Get the state trajectory
    Function xfcn = new Function("xfcn", List.of(u), List.of(s_all, v_all, m_all));
    List<DM> xopt = xfcn.call(List.of(new DM(uopt)));
    double[] sopt = xopt.get(0).nonzeros(), vopt = xopt.get(1).nonzeros(), mopt = xopt.get(2).nonzeros();
    System.out.println("position: " + str(sopt));
    System.out.println("velocity: " + str(vopt));
    System.out.println("mass:     " + str(mopt));

    // Create Matlab script to plot the solution
    String filename = "rocket_ipopt_results.m";
    PrintWriter file = new PrintWriter(filename);
    Date now = new Date();
    file.println("% Results file from rocket_ipopt.java");
    file.println("% Generated " + new SimpleDateFormat("MMM dd yyyy", Locale.US).format(now)
        + " at " + new SimpleDateFormat("HH:mm:ss").format(now));
    file.println();
    file.println("cost = " + str(cost) + ";");
    file.println("u = " + str(uopt) + ";");

    // Save results to file
    file.println("t = linspace(0,10.0," + nu + ");");
    file.println("s = " + str(sopt) + ";");
    file.println("v = " + str(vopt) + ";");
    file.println("m = " + str(mopt) + ";");

    // Finalize the results file
    file.println();
    file.println("% Plot the results");
    file.println("figure(1);");
    file.println("clf;");
    file.println();

    file.println("subplot(2,2,1);");
    file.println("plot(t,s);");
    file.println("grid on;");
    file.println("xlabel('time [s]');");
    file.println("ylabel('position [m]');");
    file.println();

    file.println("subplot(2,2,2);");
    file.println("plot(t,v);");
    file.println("grid on;");
    file.println("xlabel('time [s]');");
    file.println("ylabel('velocity [m/s]');");
    file.println();

    file.println("subplot(2,2,3);");
    file.println("plot(t,m);");
    file.println("grid on;");
    file.println("xlabel('time [s]');");
    file.println("ylabel('mass [kg]');");
    file.println();

    file.println("subplot(2,2,4);");
    file.println("plot(t,u);");
    file.println("grid on;");
    file.println("xlabel('time [s]');");
    file.println("ylabel('Thrust [kg m/s^2]');");
    file.println();

    file.close();
    System.out.println("Results saved to \"" + filename + "\"");
  }

  // C++ str(): 6 significant digits, as std::ostream
  static String str(double x) { return new DM(x).toString(); }
  static String str(double[] x) {
    StringBuilder r = new StringBuilder("[");
    for (int i = 0; i < x.length; ++i) r.append(i > 0 ? ", " : "").append(str(x[i]));
    return r.append("]").toString();
  }
}
