package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

/** Port of test/python/simulator.py */
public class SimulatorTest extends CasadiTestCase {
  private Map<String, SX> dae;
  private Function integrator, qe;
  private MX qend, q0, par;
  private final double numTend = 2.3, numQ0 = 7.1, numP = 2;

  private static double[] linspace(double a, double b, int n) {
    double[] r = new double[n];
    for (int i = 0; i < n; i++) r[i] = a + (b - a) * i / (n - 1);
    return r;
  }

  private static Map<String, Object> opts() {
    Map<String, Object> opts = new LinkedHashMap<>();
    opts.put("reltol", 1e-15);
    opts.put("abstol", 1e-15);
    opts.put("fsens_err_con", true);
    return opts;
  }

  /** helpers.checkfunction with adj=False, jacobian=False, sens_der=False, evals=False: outputs and forward */
  private void checkOutputsForward(Function trial, Function solution, List<DM> inputs, double tol) {
    assertEquals(solution.n_in(), trial.n_in());
    assertEquals(solution.n_out(), trial.n_out());
    List<DM> a = trial.call(inputs), b = solution.call(inputs);
    for (int i = 0; i < a.size(); i++) {
      if (trial.sparsity_out(i).is_empty() || solution.sparsity_out(i).is_empty()) continue;
      checkarray(a.get(i), b.get(i), tol, "output " + trial.name_out(i));
    }
    List<DM> fwd = new ArrayList<>(inputs);
    fwd.addAll(a);
    for (DM s : randomInputs(trial)) fwd.add(s);
    List<DM> fa = trial.forward(1).call(fwd), fb = solution.forward(1).call(fwd);
    for (int i = 0; i < fa.size(); i++) {
      if (fa.get(i).is_empty() || fb.get(i).is_empty()) continue;
      checkarray(fa.get(i), fb.get(i), tol, "forward " + i);
    }
  }

  @BeforeEach
  public void setUp() {
    Assumptions.assumeTrue(hasIntegrator("cvodes"));
    SX t = SX.sym("t"), q = SX.sym("q"), p = SX.sym("p");
    Map<String, SX> f = new LinkedHashMap<>();
    f.put("t", t);
    f.put("x", q);
    f.put("p", p);
    f.put("ode", q.div(p).times(ca.power(t, 2.0)));
    integrator = ca.integrator("integrator", "cvodes", f, 0, 2.3, opts());
    q0 = MX.sym("q0");
    par = MX.sym("p");
    qend = integrator.call(Map.of("x0", q0, "p", par)).get("xf");
    qe = new Function("qe", List.of(q0, par), List.of(qend));
    dae = f;
  }

  @Test
  public void test_sim_full() {
    message("Simulator inputs");
    int N = 4;
    double[] tc = linspace(0, numTend, N);
    SX t = SX.sym("t"), q = SX.sym("q"), p = SX.sym("p");
    Map<String, SX> f = new LinkedHashMap<>();
    f.put("t", t);
    f.put("x", q);
    f.put("p", p);
    f.put("ode", q.div(p).times(ca.power(t, 2.0)));
    Function integrator = ca.integrator("integrator", "cvodes", f, tc[0], tc, opts());
    List<SX> cols = new ArrayList<>();
    for (double ti : tc) cols.add(q.times(ca.exp(ca.rdivide(Math.pow(ti, 3), p.times(3.0)))));
    Map<String, SX> sol = new LinkedHashMap<>();
    sol.put("x0", q);
    sol.put("p", p);
    sol.put("xf", ca.horzcat(cols));
    Function solution = new Function("solution", sol, ca.integrator_in(), ca.integrator_out());
    Map<String, DM> fin = new LinkedHashMap<>();
    fin.put("x0", new DM(0.3));
    fin.put("p", new DM(0.7));
    List<DM> inputs = new ArrayList<>();
    for (String n : integrator.name_in()) inputs.add(fin.getOrDefault(n, DM.zeros(integrator.sparsity_in(n))));
    checkOutputsForward(integrator, solution, inputs, 1e-6);
  }

  @Test
  public void test_simulator_time_offset() {
    message("CVodes integration: simulator time offset");
    double[] t = linspace(0.7, numTend, 100);
    Function integrator = ca.integrator("integrator", "cvodes", dae, t[0], t, opts());
    int x0_ind = ca.integrator_in().indexOf("x0");
    List<DM> integrator_in = new ArrayList<>();
    for (long i = 0; i < integrator.n_in(); i++) integrator_in.add(new DM(0));
    integrator_in.set(x0_ind, dm(numQ0));
    int p_ind = ca.integrator_in().indexOf("p");
    integrator_in.set(p_ind, dm(numP));
    List<DM> integrator_out = integrator.call(integrator_in);
    double tend = numTend, q0 = numQ0, p = numP;
    assertEquals(q0 * Math.exp((Math.pow(tend, 3) - Math.pow(0.7, 3)) / (3 * p)), integrator_out.get(0).get(0, -1).to_double(), 0.5e-9, "Evaluation output mismatch");
  }
}
