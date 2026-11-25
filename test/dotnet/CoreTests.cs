using System;
using System.Collections.Generic;
using System.Linq;

namespace CasADi.Tests {

// DM values, conversions from C# numbers and arrays, indexing
public class CoreMatrixTest {
  public void ScalarAndArrays() {
    DM a = 3.0;
    Check.Near(3.0, (double)a);
    DM v = new double[] { 1, 2, 3 };
    Check.Equal((3L, 1L), v.size());
    DM m = new double[,] { { 1, 2, 3 }, { 4, 5, 6 } };
    Check.Equal((2L, 3L), m.size());
    Check.Near(new double[,] { { 1, 2, 3 }, { 4, 5, 6 } }, m);
    DM empty = new double[0];
    Check.Equal((0L, 1L), empty.size());
  }

  public void Constructors() {
    Check.Equal((2L, 3L), ca.DM(2, 3).size());
    Check.Equal(0L, ca.DM(2, 3).nnz());
    var sp = Sparsity.dense(2, 2);
    Check.Equal(4L, ca.DM(sp).nnz());
    Check.Near(new double[,] { { 7, 7 }, { 7, 7 } }, ca.DM(sp, 7.0));
    Check.Near(new double[,] { { 1, 0 }, { 0, 1 } }, DM.eye(2));
  }

  public void Arithmetic() {
    DM x = new double[] { 1, 2, 3 };
    Check.Near(new double[] { 2, 4, 6 }, x + x);
    Check.Near(new double[] { 0, 1, 2 }, x - 1);
    Check.Near(new double[] { 1, 4, 9 }, x * x);
    Check.Near(new double[] { 0.5, 1, 1.5 }, x / 2);
    Check.Near(new double[] { -1, -2, -3 }, -x);
    Check.Near(new double[] { 1, 0, 0 }, x < 2);
    Check.Near(14.0, (double)ca.mtimes(x.T(), x));
    Check.Near(Math.Sin(2.0), (double)ca.sin(ca.DM(2.0)));
    Check.Near(Math.Sin(2.0), ca.sin(2.0));
  }

  public void Indexing() {
    DM m = new double[,] { { 1, 2, 3 }, { 4, 5, 6 } };
    Check.Near(5.0, (double)m[1, 1]);
    Check.Near(4.0, (double)m[1]);          // column-major, 0-based
    Check.Near(6.0, (double)m[^1]);         // from the end
    Check.Near(new double[,] { { 2, 3 } }, m[0, 1..]);
    Check.Near(new double[,] { { 1 }, { 4 } }, m[.., 0]);
    m[0, 0] = 10;
    Check.Near(10.0, (double)m[0, 0]);
    m[1, ..] = new double[,] { { 7, 8, 9 } };
    Check.Near(new double[,] { { 10, 2, 3 }, { 7, 8, 9 } }, m);
    Check.Near(new double[] { 10, 8 }, m[new long[] { 0, 3 }]);
  }

  public void Concatenation() {
    DM a = new double[] { 1, 2 };
    Check.Equal((4L, 1L), ca.vertcat(a, a).size());
    Check.Equal((2L, 2L), ca.horzcat(a, a).size());
    Check.Equal((0L, 1L), ca.vertcat().size());
    Check.Equal((4L, 1L), ca.vertcat(new List<DM> { a, a }).size());
  }

  public void ToStringAndNonzeros() {
    DM v = new double[] { 1, 2 };
    Check.Equal("[1, 2]", v.ToString());
    Check.True(v.nonzeros().SequenceEqual(new[] { 1.0, 2.0 }));
    Check.True(Sparsity.dense(2, 2).ToString().Contains("2x2"));
    // Every printable proxy shows casadi's description, not the .NET type name
    Check.True(ca.Opti().ToString().Contains("variables"), ca.Opti().ToString());
    SX px = SX.sym("x");
    Check.True(ca.Function("f", new[] { px }, new[] { 2 * px }).ToString().StartsWith("f:"));
    Check.True(!ca.DaeBuilder("d").ToString().StartsWith("CasADi."), ca.DaeBuilder("d").ToString());
  }
}

// SX/MX expressions: overloads pick the type from the C# argument types
public class CoreSymbolicTest {
  public void SXFunction() {
    SX x = SX.sym("x");
    SX y = SX.sym("y", 2);
    SX e = ca.sin(x) * y + 1;                 // double -> SX, SX * SX
    var f = ca.Function("f", new[] { x, y }, new[] { e });
    Check.Equal(2L, f.n_in());
    var r = f.call(new DM[] { 0.5, new double[] { 1, 2 } });
    Check.Near(new[] { Math.Sin(0.5) + 1, 2 * Math.Sin(0.5) + 1 }, r[0]);
    Check.Near(new[] { Math.Sin(0.5) + 1, 2 * Math.Sin(0.5) + 1 }, f.Call(0.5, new double[] { 1, 2 }));
  }

  public void MXFunctionAndJacobian() {
    MX x = MX.sym("x", 2);
    MX e = ca.mtimes(x.T(), x);
    var f = ca.Function("f", new[] { x }, new[] { e, ca.jacobian(e, x) });
    var r = f.call(new DM[] { new double[] { 1, 2 } });
    Check.Near(5.0, (double)r[0]);
    Check.Near(new double[,] { { 2, 4 } }, r[1]);
  }

  public void MixedTypesPromote() {
    SX x = SX.sym("x");
    DM d = 2.0;
    SX e = d * x;                          // DM -> SX, as Python
    Check.True(e.is_symbolic() == false);
    MX m = MX.sym("m");
    MX em = m + d;
    Check.Equal((1L, 1L), em.size());
  }

  public void DictCall() {
    SX x = SX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { 2 * x }, new[] { "x" }, new[] { "y" });
    var r = f.call(new Dictionary<string, DM> { ["x"] = 3.0 });
    Check.Near(6.0, (double)r["y"]);
  }

  public void Options() {
    SX x = SX.sym("x");
    var opts = new Dictionary<string, object> { ["verbose"] = false, ["max_io"] = 3L };
    var f = ca.Function("f", new[] { x }, new[] { x }, opts);
    Check.Equal("f", f.name());
  }
}

// Errors from C++ arrive as SwigException; a wrong argument names the prototypes
public class CoreErrorTest {
  public void CppException() {
    var e = Check.Throws<SwigException>(() => ca.mtimes(ca.DM(2, 3), ca.DM(2, 3)));
    Check.True(e.Message.Contains("incompatible dimensions"), e.Message);
  }

  public void DisposedProxy() {
    DM a = 1.0;
    a.Dispose();
    Check.Throws<ObjectDisposedException>(() => a.nnz());
  }
}

// Optimization with Opti; zero-argument getters are properties
public class CoreOptiTest {
  public void Rosenbrock() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    opti.minimize(ca.pow(1 - x, 2) + 100 * ca.pow(y - x * x, 2));
    opti.subject_to(x + y <= 10);
    opti.solver("sqpmethod", new Dictionary<string, object> {
      ["qpsol"] = "qrqp", ["print_header"] = false, ["print_iteration"] = false, ["print_time"] = false,
      ["qpsol_options"] = new Dictionary<string, object> { ["print_iter"] = false, ["print_header"] = false } });
    OptiSol sol = opti.solve();
    Check.Near(1.0, (double)sol.value(x), 1e-5);
    Check.Near(1.0, (double)sol.value(y), 1e-5);
    Check.Equal(2L, opti.nx);
    Check.Equal((1L, 1L), opti.f.size());
  }
}

// A C# subclass of Callback is called from C++ (directors)
public class CoreCallbackTest {
  class Square : Callback {
    public int Calls;
    public Square() : base() {
      construct("square", new Dictionary<string, object>());
    }
    public override long get_n_in() => 1;
    public override long get_n_out() => 1;
    public override List<DM> eval(IEnumerable<DM> arg) {
      Calls++;
      DM x = arg.First();
      return new List<DM> { x * x };
    }
  }

  public void Eval() {
    var cb = new Square();
    var r = cb.call(new DM[] { 3.0 });
    Check.Near(9.0, (double)r[0]);
    Check.Equal(1, cb.Calls);
    // Inside a symbolic expression, called through C++
    MX x = MX.sym("x");
    var g = ca.Function("g", new[] { x }, new[] { cb.call(new[] { x })[0] * 2 });
    Check.Near(32.0, (double)g.Call(4.0));
    Check.Equal(2, cb.Calls);
  }

  class Thrower : Callback {
    public Thrower() { construct("thrower", new Dictionary<string, object>()); }
    public override List<DM> eval(IEnumerable<DM> arg) => throw new InvalidOperationException("boom");
  }

  public void ExceptionInOverride() {
    var cb = new Thrower();
    var e = Check.Throws<SwigException>(() => cb.call(new DM[] { 1.0 }));
    Check.True(e.Message.Contains("boom"), e.Message);
  }
}

// Proxies left to the garbage collector are deleted on an application thread
public class CoreLifetimeTest {
  public void ManyTemporaries() {
    for (int i = 0; i < 20000; ++i) {
      DM a = new double[] { i, i + 1 };
      SX s = SX.sym("s") * a;
    }
    GC.Collect();
    GC.WaitForPendingFinalizers();
    SwigRT.RunQueuedDeletes();
    DM b = new double[] { 1, 2 };
    Check.Near(new double[] { 2, 4 }, b * 2);
  }
}

}
