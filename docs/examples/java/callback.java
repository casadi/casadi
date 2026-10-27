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

//
// How to use Callback
// Joel Andersson
//

import org.casadi.*;
import java.util.*;

public class callback {

  static class MyCallback extends Callback {
    final double d;

    MyCallback(String name, double d) { this(name, d, Map.of()); }

    MyCallback(String name, double d, Map<String, ?> opts) {
      this.d = d;
      construct(name, opts);
    }

    // Number of inputs and outputs
    @Override public long get_n_in() { return 1; }
    @Override public long get_n_out() { return 1; }

    // Initialize the object
    @Override public void init() {
      System.out.println("initializing object");
    }

    // Evaluate numerically
    @Override public List<DM> eval(List<DM> arg) {
      DM x = arg.get(0);
      DM f = ca.sin(x.times(d));
      return List.of(f);
    }
  }

  // Derivates OPTION 2: Supply forward mode
  // Example from https://www.youtube.com/watch?v=mYOkLkS5yqc&t=4s

  static class Example4To3 extends Callback {
    Example4To3(String name) { this(name, Map.of()); }

    Example4To3(String name, Map<String, ?> opts) {
      construct(name, opts);
    }

    @Override public long get_n_in() { return 1; }
    @Override public long get_n_out() { return 1; }

    @Override public Sparsity get_sparsity_in(long i) {
      return Sparsity.dense(4, 1);
    }

    @Override public Sparsity get_sparsity_out(long i) {
      return Sparsity.dense(3, 1);
    }

    // Evaluate numerically
    @Override public List<DM> eval(List<DM> arg) {
      List<DM> v = ca.vertsplit(arg.get(0));
      DM a = v.get(0), b = v.get(1), c = v.get(2), d = v.get(3);
      DM ret = ca.vertcat(List.of(ca.sin(c).times(d).plus(ca.power(d, 2)), a.times(2).plus(c),
          ca.power(b, 2).plus(c.times(5))));
      return List.of(ret);
    }
  }

  static class Example4To3_Fwd extends Example4To3 {
    Callback fwd_callback;

    Example4To3_Fwd(String name) { super(name); }

    @Override public boolean has_forward(long nfwd) {
      // This example is written to work with a single forward seed vector
      // For efficiency, you may allow more seeds at once
      return nfwd == 1;
    }

    @Override public Function get_forward(long nfwd, String name, List<String> inames, List<String> onames,
        Map<String, ?> opts) {

      class ForwardFun extends Callback {
        ForwardFun() { this(Map.of()); }

        ForwardFun(Map<String, ?> opts) {
          construct(name, opts);
        }

        @Override public long get_n_in() { return 3; }
        @Override public long get_n_out() { return 1; }

        @Override public Sparsity get_sparsity_in(long i) {
          if (i == 0) { // nominal input
            return Sparsity.dense(4, 1);
          } else if (i == 1) { // nominal output
            return new Sparsity(3, 1);
          } else { // Forward seed
            return Sparsity.dense(4, 1);
          }
        }

        @Override public Sparsity get_sparsity_out(long i) {
          // Forward sensitivity
          return Sparsity.dense(3, 1);
        }

        // Evaluate numerically
        @Override public List<DM> eval(List<DM> arg) {
          List<DM> v = ca.vertsplit(arg.get(0));
          DM a = v.get(0), b = v.get(1), c = v.get(2), d = v.get(3);
          List<DM> s = ca.vertsplit(arg.get(2));
          DM a_dot = s.get(0), b_dot = s.get(1), c_dot = s.get(2), d_dot = s.get(3);
          System.out.println("Forward sweep with " + a_dot + " " + b_dot + " " + c_dot + " " + d_dot);
          DM w0 = ca.sin(c);
          DM w0_dot = ca.cos(c).times(c_dot);
          DM w1 = w0.times(d);
          DM w1_dot = w0_dot.times(d).plus(w0.times(d_dot));
          DM w2 = ca.power(d, 2);
          DM w2_dot = d_dot.times(2).times(d);
          DM r0 = w1.plus(w2);
          DM r0_dot = w1_dot.plus(w2_dot);
          DM w3 = a.times(2);
          DM w3_dot = a_dot.times(2);
          DM r1 = w3.plus(c);
          DM r1_dot = w3_dot.plus(c_dot);
          DM w4 = ca.power(b, 2);
          DM w4_dot = b_dot.times(2).times(b);
          DM w5 = w0.times(5);
          DM w5_dot = w0_dot.times(5);
          DM r2 = w4.plus(w5);
          DM r2_dot = w4_dot.plus(w5_dot);
          DM ret = ca.vertcat(List.of(r0_dot, r1_dot, r2_dot));
          return List.of(ret);
        }
      }
      // You are required to keep a reference alive to the returned Callback object
      fwd_callback = new ForwardFun();
      return fwd_callback;
    }
  }

  // Derivates OPTION 3: Supply reverse mode

  static class Example4To3_Rev extends Example4To3 {
    Callback rev_callback;

    Example4To3_Rev(String name) { super(name); }

    @Override public boolean has_reverse(long nadj) {
      // This example is written to work with a single forward seed vector
      // For efficiency, you may allow more seeds at once
      return nadj == 1;
    }

    @Override public Function get_reverse(long nadj, String name, List<String> inames, List<String> onames,
        Map<String, ?> opts) {

      class ReverseFun extends Callback {
        ReverseFun() { this(Map.of()); }

        ReverseFun(Map<String, ?> opts) {
          construct(name, opts);
        }

        @Override public long get_n_in() { return 3; }
        @Override public long get_n_out() { return 1; }

        @Override public Sparsity get_sparsity_in(long i) {
          if (i == 0) { // nominal input
            return Sparsity.dense(4, 1);
          } else if (i == 1) { // nominal output
            return new Sparsity(3, 1);
          } else { // Reverse seed
            return Sparsity.dense(3, 1);
          }
        }

        @Override public Sparsity get_sparsity_out(long i) {
          // Reverse sensitivity
          return Sparsity.dense(4, 1);
        }

        // Evaluate numerically
        @Override public List<DM> eval(List<DM> arg) {
          List<DM> v = ca.vertsplit(arg.get(0));
          DM a = v.get(0), b = v.get(1), c = v.get(2), d = v.get(3);
          List<DM> s = ca.vertsplit(arg.get(2));
          DM r0_bar = s.get(0), r1_bar = s.get(1), r2_bar = s.get(2);
          System.out.println("Reverse sweep with " + r0_bar + " " + r1_bar + " " + r2_bar);
          DM w0 = ca.sin(c);
          DM w1 = w0.times(d);
          DM w2 = ca.power(d, 2);
          DM r0 = w1.plus(w2);
          DM w3 = a.times(2);
          DM r1 = w3.plus(c);
          DM w4 = ca.power(b, 2);
          DM w5 = w0.times(5);
          DM r2 = w4.plus(w5);
          DM w4_bar = r2_bar;
          DM w5_bar = r2_bar;
          DM w0_bar = w5_bar.times(5);
          DM b_bar = b.times(2).times(w4_bar);
          DM w3_bar = r1_bar;
          DM c_bar = r1_bar;
          DM a_bar = w3_bar.times(2);
          DM w1_bar = r0_bar;
          DM w2_bar = r0_bar;
          DM d_bar = d.times(2).times(w2_bar);
          w0_bar = w0_bar.plus(w1_bar.times(d));
          d_bar = d_bar.plus(w0.times(w1_bar));
          c_bar = c_bar.plus(ca.cos(c).times(w0_bar));
          DM ret = ca.vertcat(List.of(a_bar, b_bar, c_bar, d_bar));
          return List.of(ret);
        }
      }
      // You are required to keep a reference alive to the returned Callback object
      rev_callback = new ReverseFun();
      return rev_callback;
    }
  }

  // Derivates OPTION 4: Supply full Jacobian

  static class Example4To3_Jac extends Example4To3 {
    Callback jac_callback;

    Example4To3_Jac(String name) { super(name); }

    @Override public boolean has_jacobian() { return true; }

    @Override public Function get_jacobian(String name, List<String> inames, List<String> onames,
        Map<String, ?> opts) {
      class JacFun extends Callback {
        JacFun() { this(Map.of()); }

        JacFun(Map<String, ?> opts) {
          construct(name, opts);
        }

        @Override public long get_n_in() { return 2; }
        @Override public long get_n_out() { return 1; }

        @Override public Sparsity get_sparsity_in(long i) {
          if (i == 0) { // nominal input
            return Sparsity.dense(4, 1);
          } else { // i==1: nominal output
            return new Sparsity(3, 1);
          }
        }

        @Override public Sparsity get_sparsity_out(long i) {
          return ca.sparsify(new DM(new double[][] {{0, 0, 1, 1}, {1, 0, 1, 0}, {0, 1, 1, 0}})).sparsity();
        }

        // Evaluate numerically
        @Override public List<DM> eval(List<DM> arg) {
          List<DM> v = ca.vertsplit(arg.get(0));
          DM b = v.get(1), c = v.get(2), d = v.get(3);
          DM ret = new DM(3, 4);
          ret.set(0, 2, d.times(ca.cos(c)));
          ret.set(0, 3, ca.sin(c).plus(d.times(2)));
          ret.set(1, 0, 2);
          ret.set(1, 2, 1);
          ret.set(2, 1, b.times(2));
          ret.set(2, 2, 5);
          return List.of(ret);
        }
      }

      // You are required to keep a reference alive to the returned Callback object
      jac_callback = new JacFun();
      return jac_callback;
    }
  }

  public static void main(String[] args) {
    // Use the function
    Callback f = new MyCallback("f", 0.5);
    Object res = f.invoke(2);
    System.out.println(res);

    // You may call the Callback symbolically
    MX x = MX.sym("x");
    System.out.println(f.invoke(x));

    // Derivates OPTION 1: finite-differences
    double eps = 1e-5;
    System.out.println(((DM) f.invoke(2 + eps)).minus((DM) f.invoke(2)).div(eps));

    f = new MyCallback("f", 0.5, Map.of("enable_fd", true));
    Function J = new Function("J", List.of(x), List.of(ca.jacobian((MX) f.invoke(x), x)));
    System.out.println(J.invoke(2));

    f = new Example4To3_Fwd("f");
    x = MX.sym("x", 4);
    J = new Function("J", List.of(x), List.of(ca.jacobian((MX) f.invoke(x), x)));
    System.out.println(J.invoke(new DM(new double[] {1, 2, 0, 3})));

    f = new Example4To3_Rev("f");
    x = MX.sym("x", 4);
    J = new Function("J", List.of(x), List.of(ca.jacobian((MX) f.invoke(x), x)));
    System.out.println(J.invoke(new DM(new double[] {1, 2, 0, 3})));

    f = new Example4To3_Jac("f");
    x = MX.sym("x", 4);
    J = new Function("J", List.of(x), List.of(ca.jacobian((MX) f.invoke(x), x)));
    System.out.println(J.invoke(new DM(new double[] {1, 2, 0, 3})));
  }
}
