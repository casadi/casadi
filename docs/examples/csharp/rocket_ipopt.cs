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
using System.Globalization;
using System.IO;
using CasADi;

public static class rocket_ipopt {
  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }
    Console.WriteLine("program started");

    // Dimensions
    int nu = 20;  // Number of control segments
    int nj = 100; // Number of integration steps per control segment

    // optimization variable
    SX u = SX.sym("u", nu); // control

    SX s_0 = 0; // initial position
    SX v_0 = 0; // initial speed
    SX m_0 = 1; // initial mass

    SX dt = 10.0 / (nj * nu); // time step
    SX alpha = 0.05; // friction
    SX beta = 0.1; // fuel consumption rate

    // Trajectory
    List<SX> s_k = new(), v_k = new(), m_k = new();

    // Integrate over the interval with Euler forward
    SX s = s_0, v = v_0, m = m_0;
    for (int k = 0; k < nu; ++k) {
      for (int j = 0; j < nj; ++j) {
        s += dt * v;
        v += dt / m * (u[k] - alpha * v * v);
        m += -dt * beta * u[k] * u[k];
      }
      s_k.Add(s);
      v_k.Add(v);
      m_k.Add(m);
    }
    SX s_all = ca.vertcat(s_k), v_all = ca.vertcat(v_k), m_all = ca.vertcat(m_k);

    // Objective function
    SX f = ca.dot(u, u);

    // Terminal constraints
    SX g = ca.vertcat(s, v, v_all);

    // Create the NLP
    var nlp = new Dictionary<string, SX> { ["x"] = u, ["f"] = f, ["g"] = g };

    // Allocate an NLP solver and buffers
    Function solver = ca.nlpsol("solver", "ipopt", nlp);

    // Bounds on g
    var gmin = new double[2 + nu];
    var gmax = new double[2 + nu];
    gmin[0] = gmax[0] = 10;
    Array.Fill(gmin, double.NegativeInfinity, 2, nu);
    Array.Fill(gmax, 1.1, 2, nu);

    // Solve the problem
    var res = solver.call(new Dictionary<string, DM> {
      ["lbx"] = -10,
      ["ubx"] = 10,
      ["x0"] = 0.4,
      ["lbg"] = gmin,
      ["ubg"] = gmax });

    // Print the optimal cost
    DM cost = res["f"];
    Console.WriteLine($"optimal cost: {cost}");

    // Print the optimal solution
    DM uopt = res["x"];
    Console.WriteLine($"optimal control: {uopt}");

    // Get the state trajectory
    var xfcn = ca.Function("xfcn", [u], [s_all, v_all, m_all]);
    List<DM> xopt = xfcn.call(new[] { uopt });
    DM sopt = xopt[0], vopt = xopt[1], mopt = xopt[2];
    Console.WriteLine($"position: {sopt}");
    Console.WriteLine($"velocity: {vopt}");
    Console.WriteLine($"mass:     {mopt}");

    // Create Matlab script to plot the solution
    string filename = "rocket_ipopt_results.m";
    using (var file = new StreamWriter(filename)) {
      DateTime now = DateTime.Now;
      file.WriteLine("% Results file from rocket_ipopt.cs");
      file.WriteLine("% Generated " + now.ToString("MMM dd yyyy", CultureInfo.InvariantCulture)
          + " at " + now.ToString("HH:mm:ss", CultureInfo.InvariantCulture));
      file.WriteLine();
      file.WriteLine($"cost = {cost};");
      file.WriteLine($"u = {uopt};");

      // Save results to file
      file.WriteLine($"t = linspace(0,10.0,{nu});");
      file.WriteLine($"s = {sopt};");
      file.WriteLine($"v = {vopt};");
      file.WriteLine($"m = {mopt};");

      // Finalize the results file
      file.WriteLine();
      file.WriteLine("% Plot the results");
      file.WriteLine("figure(1);");
      file.WriteLine("clf;");
      file.WriteLine();

      file.WriteLine("subplot(2,2,1);");
      file.WriteLine("plot(t,s);");
      file.WriteLine("grid on;");
      file.WriteLine("xlabel('time [s]');");
      file.WriteLine("ylabel('position [m]');");
      file.WriteLine();

      file.WriteLine("subplot(2,2,2);");
      file.WriteLine("plot(t,v);");
      file.WriteLine("grid on;");
      file.WriteLine("xlabel('time [s]');");
      file.WriteLine("ylabel('velocity [m/s]');");
      file.WriteLine();

      file.WriteLine("subplot(2,2,3);");
      file.WriteLine("plot(t,m);");
      file.WriteLine("grid on;");
      file.WriteLine("xlabel('time [s]');");
      file.WriteLine("ylabel('mass [kg]');");
      file.WriteLine();

      file.WriteLine("subplot(2,2,4);");
      file.WriteLine("plot(t,u);");
      file.WriteLine("grid on;");
      file.WriteLine("xlabel('time [s]');");
      file.WriteLine("ylabel('Thrust [kg m/s^2]');");
      file.WriteLine();
    }
    Console.WriteLine($"Results saved to \"{filename}\"");
  }
}
