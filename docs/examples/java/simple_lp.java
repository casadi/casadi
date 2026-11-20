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
import java.util.Map;

// minimize    3x + 4y
// subject to  x + 2y <= 14
//            3x -  y >= 0
//             x -  y <= 2

public class simple_lp {
  public static void main(String[] args) {
    if (!ca.has_conic("qpoases")) { System.out.println("SKIP: conic qpoases not available"); return; }
    double inf = Double.POSITIVE_INFINITY;

    // Sparsity of the LP linear term
    Sparsity A = Sparsity.dense(3, 2);

    // Create solver
    Function solver = ca.conic("solver", "qpoases", Map.of("a", A));
    //Function solver = ca.conic("solver", "clp", Map.of("a", A)); // Use clp

    DM g = new DM(new double[] {3, 4});
    DM a = new DM(new double[][] {{1, 2}, {3, -1}, {1, -1}});
    DM lba = new DM(new double[] {-inf, 0, -inf});
    DM uba = new DM(new double[] {14, inf, 2});

    Map<String, DM> sol = solver.call(Map.of("g", g, "a", a, "lba", lba, "uba", uba));
    System.out.println(sol);
  }
}
