package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BinaryOperator;
import java.util.function.UnaryOperator;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

/** Port of test/python/symbolic.py */
public class SymbolicTest extends CasadiTestCase {

  /** An expression in a and b, evaluable on SX, MX and double (Python builds strings and eval's them) */
  private interface Expr {
    SX sx(SX a, SX b);

    MX mx(MX a, MX b);

    double num(double a, double b);

    String str();
  }

  private static Expr var(int i) {
    return new Expr() {
      public SX sx(SX a, SX b) { return i == 0 ? a : b; }
      public MX mx(MX a, MX b) { return i == 0 ? a : b; }
      public double num(double a, double b) { return i == 0 ? a : b; }
      public String str() { return i == 0 ? "a" : "b"; }
    };
  }

  private static Expr neg(Expr x) {
    return new Expr() {
      public SX sx(SX a, SX b) { return x.sx(a, b).unaryMinus(); }
      public MX mx(MX a, MX b) { return x.mx(a, b).unaryMinus(); }
      public double num(double a, double b) { return -x.num(a, b); }
      public String str() { return "(-" + x.str() + ")"; }
    };
  }

  private static Expr bin(char op, Expr x, Expr y) {
    return new Expr() {
      public SX sx(SX a, SX b) {
        SX u = x.sx(a, b), v = y.sx(a, b);
        return op == '+' ? u.plus(v) : op == '-' ? u.minus(v) : op == '/' ? u.div(v) : u.times(v);
      }
      public MX mx(MX a, MX b) {
        MX u = x.mx(a, b), v = y.mx(a, b);
        return op == '+' ? u.plus(v) : op == '-' ? u.minus(v) : op == '/' ? u.div(v) : u.times(v);
      }
      public double num(double a, double b) {
        double u = x.num(a, b), v = y.num(a, b);
        if (op == '/' && v == 0) throw new ArithmeticException("ZeroDivisionError");
        return op == '+' ? u + v : op == '-' ? u - v : op == '/' ? u / v : u * v;
      }
      public String str() { return "(" + x.str() + op + y.str() + ")"; }
    };
  }

  @Test
  public void test_in_place_simplification() {
    Assumptions.assumeTrue(RUN_SLOW);
    System.out.println("In place simplification");

    double[] numbers = {1.3, 0.7001};
    List<BinaryOperator<Expr>> binary_operations = List.of((x, y) -> bin('+', x, y), (x, y) -> bin('-', x, y),
        (x, y) -> bin('/', x, y), (x, y) -> bin('*', x, y));
    List<UnaryOperator<Expr>> unary_operations = List.of(SymbolicTest::neg, x -> x);

    List<Expr> nodes = new ArrayList<>();
    for (int b = 0; b < 2; b++)
      for (UnaryOperator<Expr> u : unary_operations) nodes.add(u.apply(var(b)));

    List<Expr> operations_node = new ArrayList<>();
    for (Expr b1 : nodes)
      for (Expr b2 : nodes) {
        for (BinaryOperator<Expr> op : binary_operations)
          for (UnaryOperator<Expr> u : unary_operations) operations_node.add(u.apply(op.apply(b1, b2)));
        for (UnaryOperator<Expr> u : unary_operations) {
          operations_node.add(u.apply(b1));
          operations_node.add(u.apply(b2));
        }
      }
    operations_node.addAll(nodes);

    for (int t = 0; t < 2; t++) {
      SX as = SX.sym("a"), bs = SX.sym("b");
      MX am = MX.sym("a"), bm = MX.sym("b");
      int i = 0;
      for (Expr op1 : operations_node)
        for (Expr op2 : operations_node)
          for (BinaryOperator<Expr> top : binary_operations) {
            Expr op3 = top.apply(op1, op2);
            i += 1;
            Function f = t == 0 ? new Function("f", List.of(as, bs), List.of(op3.sx(as, bs)))
                : new Function("f", List.of(am, bm), List.of(op3.mx(am, bm)));
            double r = ((DM) f.invoke(numbers[0], numbers[1])).to_double();
            Double num;
            try {
              num = op3.num(numbers[0], numbers[1]);
            } catch (ArithmeticException e) {
              num = null;
            }
            if (num == null) {
              assertTrue(Double.isNaN(r) || Double.isInfinite(r));
              continue;
            }
            if (Math.abs(num - r) > 1e-10) {
              System.out.println(i + " " + op3.str() + " -> " + (t == 0 ? op3.sx(as, bs) : op3.mx(am, bm)));
              fail();
            }
          }
    }
  }
}
