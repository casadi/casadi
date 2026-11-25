using System;
using System.Collections.Generic;

namespace CasADi.Tests {

/// <summary>Port of test/java SymbolicTest (test/python/symbolic.py)</summary>
public class SymbolicTest : CasadiTestCase {
  /// <summary>An expression in a and b, evaluable on SX, MX and double</summary>
  interface Expr {
    SX sx(SX a, SX b);
    MX mx(MX a, MX b);
    double num(double a, double b);
    string str();
  }

  class Var : Expr {
    readonly int i;
    public Var(int i) { this.i = i; }
    public SX sx(SX a, SX b) => i == 0 ? a : b;
    public MX mx(MX a, MX b) => i == 0 ? a : b;
    public double num(double a, double b) => i == 0 ? a : b;
    public string str() => i == 0 ? "a" : "b";
  }

  class Neg : Expr {
    readonly Expr x;
    public Neg(Expr x) { this.x = x; }
    public SX sx(SX a, SX b) => -x.sx(a, b);
    public MX mx(MX a, MX b) => -x.mx(a, b);
    public double num(double a, double b) => -x.num(a, b);
    public string str() => "(-" + x.str() + ")";
  }

  class Bin : Expr {
    readonly char op;
    readonly Expr x, y;
    public Bin(char op, Expr x, Expr y) { this.op = op; this.x = x; this.y = y; }
    public SX sx(SX a, SX b) {
      SX u = x.sx(a, b), v = y.sx(a, b);
      return op == '+' ? u + v : op == '-' ? u - v : op == '/' ? u / v : u * v;
    }
    public MX mx(MX a, MX b) {
      MX u = x.mx(a, b), v = y.mx(a, b);
      return op == '+' ? u + v : op == '-' ? u - v : op == '/' ? u / v : u * v;
    }
    public double num(double a, double b) {
      double u = x.num(a, b), v = y.num(a, b);
      if (op == '/' && v == 0) throw new DivideByZeroException("ZeroDivisionError");
      return op == '+' ? u + v : op == '-' ? u - v : op == '/' ? u / v : u * v;
    }
    public string str() => "(" + x.str() + op + y.str() + ")";
  }

  public void test_in_place_simplification() {
    Check.Assume(RunSlow, "slow");
    Console.WriteLine("In place simplification");

    double[] numbers = { 1.3, 0.7001 };
    var binary_operations = new List<Func<Expr, Expr, Expr>> {
      (x, y) => new Bin('+', x, y), (x, y) => new Bin('-', x, y),
      (x, y) => new Bin('/', x, y), (x, y) => new Bin('*', x, y) };
    var unary_operations = new List<Func<Expr, Expr>> { x => new Neg(x), x => x };

    var nodes = new List<Expr>();
    for (int b = 0; b < 2; b++)
      foreach (var u in unary_operations) nodes.Add(u(new Var(b)));

    var operations_node = new List<Expr>();
    foreach (Expr b1 in nodes)
      foreach (Expr b2 in nodes) {
        foreach (var op in binary_operations)
          foreach (var u in unary_operations) operations_node.Add(u(op(b1, b2)));
        foreach (var u in unary_operations) {
          operations_node.Add(u(b1));
          operations_node.Add(u(b2));
        }
      }
    operations_node.AddRange(nodes);

    for (int t = 0; t < 2; t++) {
      SX as_ = SX.sym("a"), bs = SX.sym("b");
      MX am = MX.sym("a"), bm = MX.sym("b");
      int i = 0;
      foreach (Expr op1 in operations_node)
        foreach (Expr op2 in operations_node)
          foreach (var top in binary_operations) {
            Expr op3 = top(op1, op2);
            i += 1;
            Function f = t == 0 ? ca.Function("f", new[] { as_, bs }, new[] { op3.sx(as_, bs) })
                : ca.Function("f", new[] { am, bm }, new[] { op3.mx(am, bm) });
            double r = (double)f.Call(numbers[0], numbers[1]);
            double? num;
            try {
              num = op3.num(numbers[0], numbers[1]);
            } catch (DivideByZeroException) {
              num = null;
            }
            if (num == null) {
              assertTrue(double.IsNaN(r) || double.IsInfinity(r));
              continue;
            }
            if (Math.Abs(num.Value - r) > 1e-10) {
              string e = t == 0 ? op3.sx(as_, bs).ToString() : op3.mx(am, bm).ToString();
              throw new AssertFailed(i + " " + op3.str() + " -> " + e);
            }
          }
    }
  }
}

}
