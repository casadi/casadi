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

public static class rocket_snopt {
  public static void Main() {
    if (Environment.GetEnvironmentVariable("SKIP_SNOPT_TESTS") != null) return;

    if (!ca.has_nlpsol("snopt")) { Console.WriteLine("SKIP: nlpsol snopt not available"); return; }
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
    SX s_traj = SX.zeros(nu);
    SX v_traj = SX.zeros(nu);
    SX m_traj = SX.zeros(nu);

    // Integrate over the interval with Euler forward
    SX s = s_0, v = v_0, m = m_0;
    for (int k = 0; k < nu; ++k) {
      for (int j = 0; j < nj; ++j) {
        s += dt * v;
        v += dt / m * (u[k] - alpha * v * v);
        m += -dt * beta * u[k] * u[k];
      }
      s_traj[k] = s;
      v_traj[k] = v;
      m_traj[k] = m;
    }

    // Objective function
    SX f = ca.dot(u, u);

    // Terminal constraints
    SX g = ca.vertcat(s, v, v_traj);

    // Create an NLP
    var nlp = new Dictionary<string, SX> { ["x"] = u, ["f"] = f, ["g"] = g };

    // Create an NLP solver and buffers
    Function solver = ca.nlpsol("solver", "snopt", nlp);
    var arg = new Dictionary<string, DM>();

    // Bounds on u and initial condition
    arg["lbx"] = -10;
    arg["ubx"] = 10;
    arg["x0"] = 0.4;

    // Bounds on g
    var gmin = new double[2 + nu];
    var gmax = new double[2 + nu];
    gmin[0] = gmax[0] = 10;
    gmin[1] = gmax[1] = 0;
    Array.Fill(gmin, double.NegativeInfinity, 2, nu);
    Array.Fill(gmax, 1.1, 2, nu);
    arg["lbg"] = gmin;
    arg["ubg"] = gmax;

    // Solve the problem
    var res = solver.call(arg);

    // Print the optimal cost
    DM cost = res["f"];
    Console.WriteLine($"optimal cost: {cost}");

    // Print the optimal solution
    DM uopt = res["x"];
    Console.WriteLine($"optimal control: {uopt}");

    // Get the state trajectory
    var xfcn = ca.Function("xfcn", [u], [s_traj, v_traj, m_traj]);
    List<DM> xopt = xfcn.call(new[] { uopt });
    DM sopt = xopt[0], vopt = xopt[1], mopt = xopt[2];
    Console.WriteLine($"position: {sopt}");
    Console.WriteLine($"velocity: {vopt}");
    Console.WriteLine($"mass:     {mopt}");

    // Create Matlab script to plot the solution
    string filename = "rocket_snopt_results.m";
    using (var file = new StreamWriter(filename)) {
      DateTime now = DateTime.Now;
      file.WriteLine("% Results file from rocket_snopt.cs");
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
