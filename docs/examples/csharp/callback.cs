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
using System.Linq;
using CasADi;

public static class callback {

  class MyCallback : Callback {
    readonly double d;

    public MyCallback(string name, double d) : this(name, d, new Dictionary<string, object>()) { }

    public MyCallback(string name, double d, IDictionary<string, object> opts) {
      this.d = d;
      construct(name, opts);
    }

    // Number of inputs and outputs
    public override long get_n_in() => 1;
    public override long get_n_out() => 1;

    // Initialize the object
    public override void init() {
      Console.WriteLine("initializing object");
    }

    // Evaluate numerically
    public override List<DM> eval(IEnumerable<DM> arg) {
      DM x = arg.First();
      DM f = ca.sin(d * x);
      return new List<DM> { f };
    }
  }

  // Derivates OPTION 2: Supply forward mode
  // Example from https://www.youtube.com/watch?v=mYOkLkS5yqc&t=4s

  class Example4To3 : Callback {
    public Example4To3(string name) : this(name, new Dictionary<string, object>()) { }

    public Example4To3(string name, IDictionary<string, object> opts) {
      construct(name, opts);
    }

    public override long get_n_in() => 1;
    public override long get_n_out() => 1;

    public override Sparsity get_sparsity_in(long i) {
      return Sparsity.dense(4, 1);
    }

    public override Sparsity get_sparsity_out(long i) {
      return Sparsity.dense(3, 1);
    }

    // Evaluate numerically
    public override List<DM> eval(IEnumerable<DM> arg) {
      List<DM> v = ca.vertsplit(arg.First());
      DM a = v[0], b = v[1], c = v[2], d = v[3];
      DM ret = ca.vertcat(ca.sin(c) * d + ca.pow(d, 2), 2 * a + c, ca.pow(b, 2) + 5 * c);
      return new List<DM> { ret };
    }
  }

  class Example4To3_Fwd : Example4To3 {
    Callback fwd_callback;

    public Example4To3_Fwd(string name) : base(name) { }

    public override bool has_forward(long nfwd) {
      // This example is written to work with a single forward seed vector
      // For efficiency, you may allow more seeds at once
      return nfwd == 1;
    }

    class ForwardFun : Callback {
      public ForwardFun(string name) : this(name, new Dictionary<string, object>()) { }

      public ForwardFun(string name, IDictionary<string, object> opts) {
        construct(name, opts);
      }

      public override long get_n_in() => 3;
      public override long get_n_out() => 1;

      public override Sparsity get_sparsity_in(long i) {
        if (i == 0) { // nominal input
          return Sparsity.dense(4, 1);
        } else if (i == 1) { // nominal output
          return ca.Sparsity(3, 1);
        } else { // Forward seed
          return Sparsity.dense(4, 1);
        }
      }

      public override Sparsity get_sparsity_out(long i) {
        // Forward sensitivity
        return Sparsity.dense(3, 1);
      }

      // Evaluate numerically
      public override List<DM> eval(IEnumerable<DM> arg) {
        var args = arg.ToList();
        List<DM> v = ca.vertsplit(args[0]);
        DM a = v[0], b = v[1], c = v[2], d = v[3];
        List<DM> s = ca.vertsplit(args[2]);
        DM a_dot = s[0], b_dot = s[1], c_dot = s[2], d_dot = s[3];
        Console.WriteLine($"Forward sweep with {a_dot} {b_dot} {c_dot} {d_dot}");
        DM w0 = ca.sin(c);
        DM w0_dot = ca.cos(c) * c_dot;
        DM w1 = w0 * d;
        DM w1_dot = w0_dot * d + w0 * d_dot;
        DM w2 = ca.pow(d, 2);
        DM w2_dot = 2 * d_dot * d;
        DM r0 = w1 + w2;
        DM r0_dot = w1_dot + w2_dot;
        DM w3 = 2 * a;
        DM w3_dot = 2 * a_dot;
        DM r1 = w3 + c;
        DM r1_dot = w3_dot + c_dot;
        DM w4 = ca.pow(b, 2);
        DM w4_dot = 2 * b_dot * b;
        DM w5 = 5 * w0;
        DM w5_dot = 5 * w0_dot;
        DM r2 = w4 + w5;
        DM r2_dot = w4_dot + w5_dot;
        DM ret = ca.vertcat(r0_dot, r1_dot, r2_dot);
        return new List<DM> { ret };
      }
    }

    public override Function get_forward(long nfwd, string name, IEnumerable<string> inames, IEnumerable<string> onames,
        IDictionary<string, object> opts) {
      // You are required to keep a reference alive to the returned Callback object
      fwd_callback = new ForwardFun(name);
      return fwd_callback;
    }
  }

  // Derivates OPTION 3: Supply reverse mode

  class Example4To3_Rev : Example4To3 {
    Callback rev_callback;

    public Example4To3_Rev(string name) : base(name) { }

    public override bool has_reverse(long nadj) {
      // This example is written to work with a single forward seed vector
      // For efficiency, you may allow more seeds at once
      return nadj == 1;
    }

    class ReverseFun : Callback {
      public ReverseFun(string name) : this(name, new Dictionary<string, object>()) { }

      public ReverseFun(string name, IDictionary<string, object> opts) {
        construct(name, opts);
      }

      public override long get_n_in() => 3;
      public override long get_n_out() => 1;

      public override Sparsity get_sparsity_in(long i) {
        if (i == 0) { // nominal input
          return Sparsity.dense(4, 1);
        } else if (i == 1) { // nominal output
          return ca.Sparsity(3, 1);
        } else { // Reverse seed
          return Sparsity.dense(3, 1);
        }
      }

      public override Sparsity get_sparsity_out(long i) {
        // Reverse sensitivity
        return Sparsity.dense(4, 1);
      }

      // Evaluate numerically
      public override List<DM> eval(IEnumerable<DM> arg) {
        var args = arg.ToList();
        List<DM> v = ca.vertsplit(args[0]);
        DM a = v[0], b = v[1], c = v[2], d = v[3];
        List<DM> s = ca.vertsplit(args[2]);
        DM r0_bar = s[0], r1_bar = s[1], r2_bar = s[2];
        Console.WriteLine($"Reverse sweep with {r0_bar} {r1_bar} {r2_bar}");
        DM w0 = ca.sin(c);
        DM w1 = w0 * d;
        DM w2 = ca.pow(d, 2);
        DM r0 = w1 + w2;
        DM w3 = 2 * a;
        DM r1 = w3 + c;
        DM w4 = ca.pow(b, 2);
        DM w5 = 5 * w0;
        DM r2 = w4 + w5;
        DM w4_bar = r2_bar;
        DM w5_bar = r2_bar;
        DM w0_bar = 5 * w5_bar;
        DM b_bar = 2 * b * w4_bar;
        DM w3_bar = r1_bar;
        DM c_bar = r1_bar;
        DM a_bar = 2 * w3_bar;
        DM w1_bar = r0_bar;
        DM w2_bar = r0_bar;
        DM d_bar = 2 * d * w2_bar;
        w0_bar += w1_bar * d;
        d_bar += w0 * w1_bar;
        c_bar += ca.cos(c) * w0_bar;
        DM ret = ca.vertcat(a_bar, b_bar, c_bar, d_bar);
        return new List<DM> { ret };
      }
    }

    public override Function get_reverse(long nadj, string name, IEnumerable<string> inames, IEnumerable<string> onames,
        IDictionary<string, object> opts) {
      // You are required to keep a reference alive to the returned Callback object
      rev_callback = new ReverseFun(name);
      return rev_callback;
    }
  }

  // Derivates OPTION 4: Supply full Jacobian

  class Example4To3_Jac : Example4To3 {
    Callback jac_callback;

    public Example4To3_Jac(string name) : base(name) { }

    public override bool has_jacobian() => true;

    class JacFun : Callback {
      public JacFun(string name) : this(name, new Dictionary<string, object>()) { }

      public JacFun(string name, IDictionary<string, object> opts) {
        construct(name, opts);
      }

      public override long get_n_in() => 2;
      public override long get_n_out() => 1;

      public override Sparsity get_sparsity_in(long i) {
        if (i == 0) { // nominal input
          return Sparsity.dense(4, 1);
        } else { // i==1: nominal output
          return ca.Sparsity(3, 1);
        }
      }

      public override Sparsity get_sparsity_out(long i) {
        return ca.sparsify(new double[,] { { 0, 0, 1, 1 }, { 1, 0, 1, 0 }, { 0, 1, 1, 0 } }).sparsity();
      }

      // Evaluate numerically
      public override List<DM> eval(IEnumerable<DM> arg) {
        List<DM> v = ca.vertsplit(arg.First());
        DM b = v[1], c = v[2], d = v[3];
        DM ret = ca.DM(3, 4);
        ret[0, 2] = d * ca.cos(c);
        ret[0, 3] = ca.sin(c) + 2 * d;
        ret[1, 0] = 2;
        ret[1, 2] = 1;
        ret[2, 1] = 2 * b;
        ret[2, 2] = 5;
        return new List<DM> { ret };
      }
    }

    public override Function get_jacobian(string name, IEnumerable<string> inames, IEnumerable<string> onames,
        IDictionary<string, object> opts) {
      // You are required to keep a reference alive to the returned Callback object
      jac_callback = new JacFun(name);
      return jac_callback;
    }
  }

  public static void Main() {
    // Use the function
    Callback f = new MyCallback("f", 0.5);
    DM res = f.Call(2.0);
    Console.WriteLine(res);

    // You may call the Callback symbolically
    MX x = MX.sym("x");
    Console.WriteLine(f.Call(x));

    // Derivates OPTION 1: finite-differences
    double eps = 1e-5;
    Console.WriteLine((f.Call(2 + eps) - f.Call(2.0)) / eps);

    f = new MyCallback("f", 0.5, new Dictionary<string, object> { ["enable_fd"] = true });
    var J = ca.Function("J", [x], [ca.jacobian(f.Call(x), x)]);
    Console.WriteLine(J.Call(2.0));

    f = new Example4To3_Fwd("f");
    x = MX.sym("x", 4);
    J = ca.Function("J", [x], [ca.jacobian(f.Call(x), x)]);
    Console.WriteLine(J.Call(new double[] { 1, 2, 0, 3 }));

    f = new Example4To3_Rev("f");
    x = MX.sym("x", 4);
    J = ca.Function("J", [x], [ca.jacobian(f.Call(x), x)]);
    Console.WriteLine(J.Call(new double[] { 1, 2, 0, 3 }));

    f = new Example4To3_Jac("f");
    x = MX.sym("x", 4);
    J = ca.Function("J", [x], [ca.jacobian(f.Call(x), x)]);
    Console.WriteLine(J.Call(new double[] { 1, 2, 0, 3 }));
  }
}
