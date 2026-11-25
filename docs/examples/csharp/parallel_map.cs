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
using System.Diagnostics;
using CasADi;

public static class parallel_map {
  public static void Main() {
    Console.WriteLine("importing casadi...");

    // number of inputs to evaluate in parallel
    int N = 300;

    // dummy input
    var dummyInput = new double[N];
    for (int i = 0; i < N; i++) dummyInput[i] = 2.0 * Math.PI * i / (N - 1);

    // make a dummy function that's moderately expensive to evaluate
    Console.WriteLine("creating dummy function....");
    SX x = SX.sym("x");
    SX y = x;
    for (int k = 0; k < 100000; k++) {
      y = ca.sin(y);
    }
    var f0 = ca.Function("f", [x], [y]);

    // evaluate it serially, the old-fasioned way
    MX X = MX.sym("x", N);
    var Ys = new List<MX>();
    for (int k = 0; k < N; k++) Ys.Add(f0.Call(X[k]));
    MX Y = ca.vertcat(Ys);
    var fNaiveParallel = ca.Function("fParallel", [X], [Y]);

    Console.WriteLine("evaluating naive parallel function...");
    var sw = Stopwatch.StartNew();
    DM outNaive = fNaiveParallel.Call(dummyInput);
    Console.WriteLine($"evaluated naive parallel function in {sw.Elapsed.TotalSeconds:F3} seconds");

    // evaluate it using new serial map construct
    Function fMap = f0.map(N);

    Console.WriteLine("evaluating serial map function...");
    sw.Restart();
    DM outMap = fMap.Call(dummyInput);
    Console.WriteLine($"evaluated serial map function in {sw.Elapsed.TotalSeconds:F3} seconds");
    // the following has different shaped outputs, so it's commented out
    //Console.WriteLine(outNaive == outMap);

    // evaluate it using new parallel map construct
    fMap = f0.map(N, "openmp");

    Console.WriteLine("evaluating parallel map function...");
    sw.Restart();
    outMap = fMap.Call(dummyInput);
    Console.WriteLine($"evaluated parallel map function in {sw.Elapsed.TotalSeconds:F3} seconds");
    // the following has different shaped outputs, so it's commented out
    //Console.WriteLine(outNaive == outMap);
  }
}
