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
import org.casadi.*;
import java.util.ArrayList;
import java.util.List;

public class parallel_map {
  public static void main(String[] args) {
    System.out.println("importing casadi...");

    // number of inputs to evaluate in parallel
    int N = 300;

    // dummy input
    double[] dummyInput = new double[N];
    for (int i = 0; i < N; i++) dummyInput[i] = 2.0 * Math.PI * i / (N - 1);

    // make a dummy function that's moderately expensive to evaluate
    System.out.println("creating dummy function....");
    SX x = SX.sym("x");
    SX y = x;
    for (int k = 0; k < 100000; k++) {
      y = ca.sin(y);
    }
    Function f0 = new Function("f", List.of(x), List.of(y));

    // evaluate it serially, the old-fasioned way
    MX X = MX.sym("x", N);
    List<MX> Ys = new ArrayList<>();
    for (int k = 0; k < N; k++) Ys.add((MX) f0.invoke(X.get(k)));
    MX Y = ca.vertcat(Ys);
    Function fNaiveParallel = new Function("fParallel", List.of(X), List.of(Y));

    System.out.println("evaluating naive parallel function...");
    long t0 = System.nanoTime();
    Object outNaive = fNaiveParallel.invoke(dummyInput);
    long t1 = System.nanoTime();
    System.out.println(String.format("evaluated naive parallel function in %.3f seconds", (t1 - t0) * 1e-9));

    // evaluate it using new serial map construct
    Function fMap = f0.map(N);

    System.out.println("evaluating serial map function...");
    t0 = System.nanoTime();
    Object outMap = fMap.invoke(dummyInput);
    t1 = System.nanoTime();
    System.out.println(String.format("evaluated serial map function in %.3f seconds", (t1 - t0) * 1e-9));
    // the following has different shaped outputs, so it's commented out
    //System.out.println(outNaive.equals(outMap));

    // evaluate it using new parallel map construct
    fMap = f0.map(N, "openmp");

    System.out.println("evaluating parallel map function...");
    t0 = System.nanoTime();
    outMap = fMap.invoke(dummyInput);
    t1 = System.nanoTime();
    System.out.println(String.format("evaluated parallel map function in %.3f seconds", (t1 - t0) * 1e-9));
    // the following has different shaped outputs, so it's commented out
    //System.out.println(outNaive.equals(outMap));
  }
}
