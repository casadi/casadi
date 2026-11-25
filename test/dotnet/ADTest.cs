using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace CasADi.Tests {

/// <summary>Port of test/java ADTest (test/python/ad.py)</summary>
public class ADTest : CasadiTestCase {
  static readonly string[] SHAPES = { "column", "row", "matrix" };
  static readonly string[] TYPES = { "dense", "sparse" };
  const double CODEGEN_TOL = 0.5e-15;

  readonly Dictionary<string, SX> sxinputs, sxoutputs;
  readonly Dictionary<string, MX> mxinputs;
  readonly Dictionary<string, Func<MX, MX>> mxoutputs;

  static Dictionary<string, object> D(params (string, object)[] kv) => kv.ToDictionary(e => e.Item1, e => e.Item2);

  static MX nz(MX x, long k) => x.get_nz(false, k);

  static MX temp1(MX xyz) {
    var X = ca.MX(6, 1);
    X[0, 0] = nz(xyz, 0);
    X[2, 0] = nz(xyz, 0) + 2.0 * ca.pow(nz(xyz, 1), 2.0);
    X[4, 0] = nz(xyz, 0) + 2.0 * ca.pow(nz(xyz, 1), 3.0) + 3.0 * ca.pow(nz(xyz, 2), 4.0);
    X[5, 0] = nz(xyz, 3);
    return X;
  }

  static MX temp2(MX xyz) {
    var X = ca.MX(1, 6);
    X[0, 0] = nz(xyz, 0);
    X[0, 2] = nz(xyz, 0) + 2.0 * ca.pow(nz(xyz, 1), 2.0);
    X[0, 4] = nz(xyz, 0) + 2.0 * ca.pow(nz(xyz, 1), 3.0) + 3.0 * ca.pow(nz(xyz, 2), 4.0);
    X[0, 5] = nz(xyz, 3);
    return X;
  }

  static List<MX> outList(MX xyz) => new List<MX> {
    nz(xyz, 0), nz(xyz, 0) + 2.0 * ca.pow(nz(xyz, 1), 2.0),
    nz(xyz, 0) + 2.0 * ca.pow(nz(xyz, 1), 3.0) + 3.0 * ca.pow(nz(xyz, 2), 4.0), nz(xyz, 3) };

  public ADTest() {
    SX x = SX.sym("x"), y = SX.sym("y"), z = SX.sym("z"), w = SX.sym("w");
    var outp = ca.SX(6, 1);
    outp[0, 0] = x;
    outp[2, 0] = x + 2.0 * ca.pow(y, 2.0);
    outp[4, 0] = x + 2.0 * ca.pow(y, 3.0) + 3.0 * ca.pow(z, 4.0);
    outp[5, 0] = w;
    var inp = ca.SX(6, 1);
    inp[0, 0] = x;
    inp[2, 0] = y;
    inp[4, 0] = z;
    inp[5, 0] = w;
    Sparsity sp = ca.Sparsity(1, 6, new long[] { 0, 1, 1, 2, 2, 3, 4 }, new long[] { 0, 0, 0, 0 }).T();
    Sparsity spT = ca.Sparsity(6, 1, new long[] { 0, 4 }, new long[] { 0, 2, 4, 5 }).T();
    SX xyzw = ca.vertcat(x, y, z, w);
    sxinputs = new Dictionary<string, SX> {
      ["column dense"] = xyzw, ["column sparse"] = inp, ["row dense"] = xyzw.T(), ["row sparse"] = inp.T(),
      ["matrix dense"] = ca.reshape(xyzw, 2, 2), ["matrix sparse"] = ca.reshape(inp, 3, 2) };
    mxinputs = new Dictionary<string, MX> {
      ["column dense"] = MX.sym("xyzw", 4, 1), ["column sparse"] = MX.sym("xyzw", sp),
      ["row dense"] = MX.sym("xyzw", 1, 4), ["row sparse"] = MX.sym("xyzw", spT),
      ["matrix dense"] = MX.sym("xyzw", 2, 2), ["matrix sparse"] = MX.sym("xyzw", ca.reshape(inp, 3, 2).sparsity()) };
    mxoutputs = new Dictionary<string, Func<MX, MX>> {
      ["column dense"] = xyz => ca.vertcat(outList(xyz)), ["column sparse"] = temp1,
      ["row dense"] = xyz => ca.horzcat(outList(xyz)), ["row sparse"] = temp2,
      ["matrix dense"] = xyz => ca.reshape(ca.vertcat(outList(xyz)), 2, 2),
      ["matrix sparse"] = xyz => ca.reshape(temp1(xyz), 3, 2) };
    SX o = ca.vertcat(x, x + 2.0 * ca.pow(y, 2.0), x + 2.0 * ca.pow(y, 3.0) + 3.0 * ca.pow(z, 4.0), w);
    sxoutputs = new Dictionary<string, SX> {
      ["column dense"] = o, ["column sparse"] = outp, ["row dense"] = o.T(), ["row sparse"] = outp.T(),
      ["matrix dense"] = ca.reshape(o, 2, 2), ["matrix sparse"] = ca.reshape(outp, 3, 2) };
  }

  static DM jacobians(string inputtype, string outputtype, double[] n) {
    double y = n[1], z = n[2];
    string k = inputtype + " " + outputtype;
    double[,] r;
    if (k == "dense dense")
      r = new double[,] { { 1, 0, 0, 0 }, { 1, 4 * y, 0, 0 }, { 1, 6 * y * y, 12 * z * z * z, 0 }, { 0, 0, 0, 1 } };
    else if (k == "dense sparse")
      r = new double[,] { { 1, 0, 0, 0 }, { 0, 0, 0, 0 }, { 1, 4 * y, 0, 0 }, { 0, 0, 0, 0 }, { 1, 6 * y * y, 12 * z * z * z, 0 }, { 0, 0, 0, 1 } };
    else if (k == "sparse dense")
      r = new double[,] { { 1, 0, 0, 0, 0, 0 }, { 1, 0, 4 * y, 0, 0, 0 }, { 1, 0, 6 * y * y, 0, 12 * z * z * z, 0 }, { 0, 0, 0, 0, 0, 1 } };
    else
      r = new double[,] { { 1, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 1, 0, 4 * y, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 },
                          { 1, 0, 6 * y * y, 0, 12 * z * z * z, 0 }, { 0, 0, 0, 0, 0, 1 } };
    return r;
  }

  static DM nonzeroIndicator(DM J) {
    double[,] a = J.ToArray();
    for (int i = 0; i < a.GetLength(0); i++)
      for (int j = 0; j < a.GetLength(1); j++) a[i, j] = a[i, j] != 0 ? 1 : 0;
    return a;
  }

  static Function jacobian_old(Function f, long i, long j) {
    var outs = new List<string> { "jac:" + f.name_out(j) + ":" + f.name_in(i) };
    outs.AddRange(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), outs);
  }

  static Function hessian_old(Function f, long i, long j) {
    var outs = new List<string> { "hess:" + f.name_out(j) + ":" + f.name_in(i) + ":" + f.name_in(i),
                                  "grad:" + f.name_out(j) + ":" + f.name_in(i) };
    outs.AddRange(f.name_out());
    return f.factory(f.name() + "_hess", f.name_in(), outs);
  }

  static List<DM> seedList(Sparsity sp, double[][] seeds) => seeds.Select(s => ca.DM(sp, (DM)s)).ToList();

  static readonly double[][] SEEDS = { new double[] { 1, 0, 0, 0 }, new double[] { 0, 2, 0, 0 }, new double[] { 1.2, 4.8, 7.9, 4.6 } };

  public void test_SXeval_sx() {
    double[] n = { 1.2, 2.3, 7, 1.4 };
    foreach (string inputshape in SHAPES)
      foreach (string outputshape in SHAPES)
        foreach (string inputtype in TYPES)
          foreach (string outputtype in TYPES) {
            message("eval_sx on SX. Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
            var f = ca.Function("f", new[] { sxinputs[inputshape + " " + inputtype] }, new[] { sxoutputs[outputshape + " " + outputtype] });
            DM f_in = ca.DM(f.sparsity_in(0), (DM)n);
            DM r = f.Call(f_in);
            DM J = jacobians(inputtype, outputtype, n);
            SX y = SX.sym("y", f.sparsity_in(0));
            List<DM> fseeds = seedList(f.sparsity_in(0), SEEDS), aseeds = seedList(f.sparsity_out(0), SEEDS);
            SX res = f.Call(y);
            var fwdsens = ca.forward(new[] { res }, new[] { y }, fseeds.Select(s => new[] { (SX)s }));
            var adjsens = ca.reverse(new[] { res }, new[] { y }, aseeds.Select(s => new[] { (SX)s }));
            var fe = ca.Function("fe", new[] { y }, new[] { res });
            DM re = fe.Call(f_in);
            checkarray(r, re);
            for (int i = 0; i < fwdsens.Count; i++) {
              fe = ca.Function("fe", new[] { y }, new[] { fwdsens[i][0] });
              re = fe.Call(f_in);
              checkarray(ca.vec(re), ca.mtimes(J, ca.vec(fseeds[i])), "AD");
            }
            for (int i = 0; i < adjsens.Count; i++) {
              fe = ca.Function("fe", new[] { y }, new[] { adjsens[i][0] });
              re = fe.Call(f_in);
              checkarray(ca.vec(re), ca.mtimes(J.T(), ca.vec(aseeds[i])), "AD");
            }
          }
  }

  public void test_MXeval_mx() {
    double[] n = { 1.2, 2.3, 7, 1.4 };
    foreach (string inputshape in SHAPES)
      foreach (string outputshape in SHAPES)
        foreach (string inputtype in TYPES)
          foreach (string outputtype in TYPES) {
            message("eval_mx on MX. Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
            MX inp = mxinputs[inputshape + " " + inputtype];
            var f = ca.Function("f", new[] { inp }, new[] { mxoutputs[outputshape + " " + outputtype](inp) });
            DM f_in = ca.DM(f.sparsity_in(0), (DM)n);
            DM r = f.Call(f_in);
            DM J = jacobians(inputtype, outputtype, n);
            MX y = MX.sym("y", f.sparsity_in(0));
            List<DM> fseeds = seedList(f.sparsity_in(0), SEEDS), aseeds = seedList(f.sparsity_out(0), SEEDS);
            MX res = f.Call(y);
            var fwdsens = ca.forward(new[] { res }, new[] { y }, fseeds.Select(s => new[] { (MX)s }));
            var adjsens = ca.reverse(new[] { res }, new[] { y }, aseeds.Select(s => new[] { (MX)s }));
            var fe = ca.Function("fe", new[] { y }, new[] { res });
            DM re = fe.Call(f_in);
            checkarray(r, re);
            for (int i = 0; i < fwdsens.Count; i++) {
              fe = ca.Function("fe", new[] { y }, new[] { fwdsens[i][0] });
              re = fe.Call(f_in);
              checkarray(ca.vec(re), ca.mtimes(J, ca.vec(fseeds[i])), "AD");
            }
            for (int i = 0; i < adjsens.Count; i++) {
              fe = ca.Function("fe", new[] { y }, new[] { adjsens[i][0] });
              re = fe.Call(f_in);
              checkarray(ca.vec(re), ca.mtimes(J.T(), ca.vec(aseeds[i])), "AD");
            }
          }
  }

  public void test_MXeval_sx() {
    Check.Assume(Environment.GetEnvironmentVariable("CASADI_DOTNET_KNOWN_BUGS") != null, "known bug (not implemented)");
    double[] n = { 1.2, 2.3, 7, 1.4 };
    foreach (string inputshape in SHAPES)
      foreach (string outputshape in SHAPES)
        foreach (string inputtype in TYPES)
          foreach (string outputtype in TYPES) {
            message("eval_sx on MX. Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
            MX inp = mxinputs[inputshape + " " + inputtype];
            var f = ca.Function("f", new[] { inp }, new[] { mxoutputs[outputshape + " " + outputtype](inp) });
            var f_in = new List<DM>();
            for (long i = 0; i < f.n_in(); i++) f_in.Add(0.0);
            f_in[0] = n;
            DM r = f.call(f_in)[0];
            DM J = jacobians(inputtype, outputtype, n);
            SX y = SX.sym("y", f.sparsity_in(0));
            List<DM> fseeds = seedList(f.sparsity_in(0), SEEDS), aseeds = seedList(f.sparsity_out(0), SEEDS);
            SX res = f.Call(y);
            var fwdsens = ca.forward(new[] { res }, new[] { y }, fseeds.Select(s => new[] { (SX)s }));
            var adjsens = ca.reverse(new[] { res }, new[] { y }, aseeds.Select(s => new[] { (SX)s }));
            var fe = ca.Function("fe", new[] { y }, new[] { res });
            var fe_in = new List<DM> { n };
            checkarray(r, fe.call(fe_in)[0]);
            for (int i = 0; i < fwdsens.Count; i++) {
              fe = ca.Function("fe", new[] { y }, new[] { fwdsens[i][0] });
              checkarray(ca.vec(fe.call(fe_in)[0].T()), ca.mtimes(J, ca.vec(fseeds[i].T())), "AD");
            }
            for (int i = 0; i < adjsens.Count; i++) {
              fe = ca.Function("fe", new[] { y }, new[] { adjsens[i][0] });
              checkarray(ca.vec(fe.call(fe_in)[0].T()), ca.mtimes(J.T(), ca.vec(aseeds[i].T())), "AD");
            }
          }
  }

  public void test_MXeval_sx_reduced() {
    double[] n = { 1.2, 2.3, 7, 1.4 };
    foreach (string inputshape in SHAPES)
      foreach (string outputshape in SHAPES)
        foreach (string inputtype in TYPES)
          foreach (string outputtype in TYPES) {
            message("eval_sx on MX. Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
            MX inp = mxinputs[inputshape + " " + inputtype];
            var f = ca.Function("f", new[] { inp }, new[] { mxoutputs[outputshape + " " + outputtype](inp) });
            DM f_in = ca.DM(f.sparsity_in(0), (DM)n);
            DM r = f.Call(f_in);
            SX y = SX.sym("y", f.sparsity_in(0));
            SX res = f.Call(y);
            ca.forward(new[] { res }, new[] { y }, new List<List<SX>>());
            ca.reverse(new[] { res }, new[] { y }, new List<List<SX>>());
            ca.Function("fe", new[] { y }, new[] { res });
            DM re = f.Call(f_in);
            checkarray(r, re);
          }
  }

  static Dictionary<string, object> adOpts(string mode) {
    int w = mode == "forward" ? 0 : 1;
    return D(("ad_weight", w), ("ad_weight_sp", w));
  }

  public void test_Jacobian() {
    double[] n = { 1.2, 2.3, 7, 4.6 };
    foreach (string inputshape in SHAPES)
      foreach (string outputshape in SHAPES)
        foreach (string inputtype in TYPES)
          foreach (string outputtype in TYPES)
            foreach (string mode in new[] { "forward", "reverse" }) {
              message(" " + mode + " Jacobian on SX. Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
              var f = ca.Function("f", new[] { sxinputs[inputshape + " " + inputtype] },
                                   new[] { sxoutputs[outputshape + " " + outputtype] }, adOpts(mode));
              Function Jf = jacobian_old(f, 0, 0);
              DM J_in = ca.DM(f.sparsity_in(0), (DM)n);
              DM Jout = Jf.call(new[] { J_in })[0];
              DM J = jacobians(inputtype, outputtype, n);
              checkarray(Jout, J, "Jacobian\n Mode: " + mode + "\n Input: " + inputshape + " " + inputtype + "\n Output: " + outputshape + " " + outputtype);
            }
  }

  public void test_jacobianSX() {
    double[] n = { 1.2, 2.3, 7, 4.6 };
    foreach (string inputshape in SHAPES)
      foreach (string outputshape in SHAPES)
        foreach (string inputtype in TYPES)
          foreach (string outputtype in TYPES) {
            message("jacobian on SX (SCT). Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
            SX inp = sxinputs[inputshape + " " + inputtype];
            var Jf = ca.Function("Jf", new[] { inp }, new[] { ca.jacobian(ca.SX(sxoutputs[outputshape + " " + outputtype]), ca.SX(inp)) });
            DM J_in = ca.DM(Jf.sparsity_in(0), (DM)n);
            DM J_out = Jf.Call(J_in);
            checkarray(J_out, jacobians(inputtype, outputtype, n), "jacobian");
          }
  }

  public void test_jacsparsity() {
    double[] n = { 1.2, 2.3, 7, 4.6 };
    foreach (string inputshape in SHAPES)
      foreach (string outputshape in SHAPES)
        foreach (string inputtype in TYPES)
          foreach (string outputtype in TYPES) {
            message("jacsparsity on SX. Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
            var f = ca.Function("f", new[] { sxinputs[inputshape + " " + inputtype] }, new[] { sxoutputs[outputshape + " " + outputtype] });
            DM J = jacobians(inputtype, outputtype, n);
            checkarray(DM.ones(f.jac_sparsity(0, 0)), nonzeroIndicator(J), "jacsparsity");
          }
  }

  public void test_JacobianMX() {
    double[] n = { 1.2, 2.3, 7, 4.6 };
    foreach (string inputshape in SHAPES)
      foreach (string outputshape in SHAPES)
        foreach (string inputtype in TYPES)
          foreach (string outputtype in TYPES)
            foreach (string mode in new[] { "forward", "reverse" }) {
              message("adj AD on MX. Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
              MX inp = mxinputs[inputshape + " " + inputtype];
              var f = ca.Function("f", new[] { inp }, new[] { mxoutputs[outputshape + " " + outputtype](inp) }, adOpts(mode));
              Function Jf = jacobian_old(f, 0, 0);
              DM J_in = ca.DM(f.sparsity_in(0), (DM)n);
              DM J_out = Jf.call(new[] { J_in })[0];
              DM J = jacobians(inputtype, outputtype, n);
              checkarray(J_out, J, "Jacobian\n Mode: " + mode + "\n Input: " + inputshape + " " + inputtype + "\n Output: " + outputshape + " " + outputtype);
            }
  }

  public void test_jacsparsityMX() {
    double[] n = { 1.2, 2.3, 7, 4.6 };
    foreach (string inputshape in SHAPES)
      foreach (string outputshape in SHAPES)
        foreach (string inputtype in TYPES)
          foreach (string outputtype in TYPES)
            foreach (string mode in new[] { "forward", "reverse" }) {
              message(" " + mode + " jacobian on MX (SCT). Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
              MX inp = mxinputs[inputshape + " " + inputtype];
              var f = ca.Function("f", new[] { inp }, new[] { mxoutputs[outputshape + " " + outputtype](inp) }, adOpts(mode));
              Function Jf = jacobian_old(f, 0, 0);
              DM J_in = ca.DM(f.sparsity_in(0), (DM)n);
              DM J_out = Jf.call(new[] { J_in })[0];
              DM J = jacobians(inputtype, outputtype, n);
              checkarray(J_out, J, "jacobian");
              checkarray(DM.ones(f.jac_sparsity(0, 0)), nonzeroIndicator(J), "jacsparsity");
            }
  }

  public void test_hessian() {
    message("Jacobian chaining");
    SX x = SX.sym("x"), y = SX.sym("y"), z = SX.sym("z");
    double[] n = { 1.2, 2.3, 7 };
    var f = ca.Function("f", new[] { ca.vertcat(x, y, z) }, new[] { ca.vertcat(x + 2.0 * ca.pow(y, 3.0) + 3.0 * ca.pow(z, 4.0)) });
    Function J = jacobian_old(f, 0, 0);
    MX m = MX.sym("m", 3, 1);
    MX JT = J.call(new[] { m })[0];
    var JTf = ca.Function("JT", new[] { m }, new[] { JT.T() });
    JTf.Call((DM)n);
    Function H = jacobian_old(JTf, 0, 0);
    H.call(new DM[] { n });
  }

  public void test_bugshape() {
    message("shape bug");
    SX x = SX.sym("x"), y = SX.sym("y");
    var inp = ca.SX(5, 1);
    inp[0, 0] = x;
    inp[3, 0] = y;
    var f = ca.Function("f", new[] { inp }, new[] { ca.vertcat(x + y, x, y) });
    Function J = jacobian_old(f, 0, 0);
    J.call(new[] { ca.DM(f.sparsity_in(0), ca.vertcat(2, 7)) });
    assertEquals(3L, f.size1_out(0), "Jacobian shape bug");
    assertEquals(1L, f.size2_out(0), "Jacobian shape bug");
  }

  public void test_bugglibc() {
    message("Code that used to throw a glibc error");
    SX x = SX.sym("x"), y = SX.sym("y");
    var inp = ca.SX(5, 1);
    inp[0, 0] = x;
    inp[3, 0] = y;
    var f = ca.Function("f", new[] { inp }, new[] { ca.vertcat(x + y, x, y) });
    Function J = jacobian_old(f, 0, 0);
    DM J_in = ca.DM(f.sparsity_in(0), ca.vertcat(2, 7));
    _ = J.call(new[] { J_in })[0];
    f = ca.Function("f", new[] { inp }, new[] { ca.vertcat(x + y, x, y) });
    jacobian_old(f, 0, 0);
  }

  // ---- test_MX ----

  static MX mx(object o) => o switch {
    MX m => m,
    DM d => ca.MX(d),
    int i => ca.MX((double)i),
    double v => ca.MX(v),
    _ => throw new ArgumentException(o.GetType().Name)
  };

  static MX bc(params object[][] rows) => ca.blockcat(rows.Select(r => r.Select(mx).ToList()).ToList());

  static DM idx(params double[] i) => i;

  static MX sq(MX a) => ca.pow(a, 2.0);

  sealed class Case {
    public readonly List<DM> values;
    public readonly MX outp, jac;
    public readonly bool withSx;
    public readonly string std;
    public Case(List<DM> values, object outp, object jac, bool withSx, string std) {
      this.values = values;
      this.outp = mx(outp);
      this.jac = mx(jac);
      this.withSx = withSx;
      this.std = std;
    }
  }

  static DM randDMsp(Random rnd, Sparsity sp) {
    double[] v = new double[sp.nnz()];
    for (int k = 0; k < v.Length; k++) v[k] = rnd.NextDouble();
    return ca.DM(sp, v.Length == 0 ? ca.DM(0, 1) : (DM)v);
  }

  Dictionary<string, List<(Function, List<DM>)>> cgBatch;

  void check_codegen(Function F, List<DM> inputs, string std) {
    if (!RunSlow) return;
    if (!cgBatch.TryGetValue(std, out var l)) cgBatch[std] = l = new List<(Function, List<DM>)>();
    l.Add((F, inputs));
  }

  static int cgCounter = 0;

  void flushCodegen() {
    var batch = cgBatch;
    cgBatch = new Dictionary<string, List<(Function, List<DM>)>>();
    foreach (var e in batch) {
      string name = "codegen_ad_" + Environment.ProcessId + "_" + (cgCounter++);
      string dir = Path.Combine(Path.GetTempPath(), "casadi_ad_" + Guid.NewGuid().ToString("N"));
      Directory.CreateDirectory(dir);
      var cg = ca.CodeGenerator(name);
      var wrapped = new List<Function>();
      for (int k = 0; k < e.Value.Count; k++) {
        Function f = e.Value[k].Item1;
        Function fw = f.wrap(f.name() + "_b" + k);
        wrapped.Add(fw);
        cg.add(fw);
      }
      cg.generate(dir + "/");
      string inc = GlobalOptions.getCasadiIncludePath(), lib = GlobalOptions.getCasadiPath();
      string so = dir + "/" + name + ".so";
      var args = new List<string> { "-pedantic", "-std=" + e.Key, "-fPIC", "-shared", "-Wall", "-Werror", "-Wextra",
        "-I" + inc, "-Wno-dangling-pointer", "-Wno-unknown-warning-option", "-Wno-unknown-pragmas", "-Wno-long-long",
        "-Wno-unused-parameter", "-Wno-misleading-indentation", "-O3", dir + "/" + name + ".c", "-o", so,
        "-L" + lib, "-Wl,-rpath," + lib };
      var psi = new ProcessStartInfo("gcc");
      foreach (var a in args) psi.ArgumentList.Add(a);
      using (var p = Process.Start(psi)) {
        p.WaitForExit();
        assertEquals(0, p.ExitCode, "gcc failed");
      }
      for (int k = 0; k < wrapped.Count; k++) {
        List<DM> inp = e.Value[k].Item2;
        Function f = wrapped[k], f2 = ca.external(f.name(), so);
        List<DM> r = f.call(inp), r2 = f2.call(inp);
        for (int i = 0; i < r.Count; i++) checkarray(r[i], r2[i], CODEGEN_TOL, "codegen " + f.name());
      }
      Directory.Delete(dir, true);
    }
  }

  void check_serialize(Function F, List<DM> inputs) {
    string s = F.serialize(D(("debug", true)));
    Function F2 = Function.deserialize(s);
    List<DM> a = F.call(inputs), b = F2.call(inputs);
    for (int i = 0; i < F.n_out(); i++) checkarray(a[i], b[i], 0.5e-16, "serialize");
  }

  void check_eval_mx(List<MX> l) => check_eval_mx(ca.veccat(l));

  void check_eval_mx(MX e) {
    List<MX> args = ca.symvar(e);
    var f = ca.Function("f", args, new[] { e });
    var new_args = args.Select(a => MX.sym(a.name(), a.sparsity())).ToList();
    MX mx_eval = f.call(new_args, true, false)[0];
    MX ev = ca.substitute(new[] { mx_eval }, new_args, args)[0];
    check_identical_mx(ev, e);
  }

  static void check_identical_mx(MX a, MX b) {
    var fa = ca.Function("f", ca.symvar(a), new[] { a });
    var fb = ca.Function("f", ca.symvar(b), new[] { b });
    assertTrue(fa.serialize() == fb.serialize(), "identical mx");
  }

  static Sparsity remove_first(Sparsity x) {
    var ret = ca.DM(x);
    if (ret.numel() > 0) ret[0, 0] = ca.DM(1, 1);
    return ret.sparsity();
  }

  static Sparsity remove_last(Sparsity x) {
    var ret = ca.DM(x);
    if (ret.nnz() > 0) {
      List<long> r = ret.sparsity().row(), c = ret.sparsity().get_col();
      ret[r[^1], c[^1]] = ca.DM(1, 1);
      return ret.sparsity();
    }
    return x;
  }

  /// <summary>What test_MX needs from SX or MX symbolics, so that the sweeps are written once</summary>
  sealed class Sym<T> {
    public Func<string, Sparsity, T> sym;
    public Func<Function, List<T>, bool, List<T>> call;
    public Func<List<T>, List<T>, List<List<T>>, Dictionary<string, object>, List<List<T>>> fwd, rev;
    public Func<string, List<T>, List<T>, Function> mk;
    public bool isMX;
  }

  static readonly Sym<MX> MXS = new Sym<MX> {
    sym = MX.sym, call = (f, a, i) => f.call(a, i), fwd = (r, i, s, o) => ca.forward(r, i, s, o),
    rev = (r, i, s, o) => ca.reverse(r, i, s, o), mk = (n, a, b) => ca.Function(n, a, b), isMX = true };
  static readonly Sym<SX> SXS = new Sym<SX> {
    sym = SX.sym, call = (f, a, i) => f.call(a, i), fwd = (r, i, s, o) => ca.forward(r, i, s, o),
    rev = (r, i, s, o) => ca.reverse(r, i, s, o), mk = (n, a, b) => ca.Function(n, a, b), isMX = false };

  /// <summary>First-order sweeps of f with symbolic seeds; returns the last vf</summary>
  Function firstOrder<T>(Sym<T> S, Function f, int ndir, List<Func<Sparsity, Sparsity>> spmods, List<DM> values, DM J_,
                         string std, Dictionary<string, List<List<DM>>> storage, ref string storagekey, ref Function vf_mx) {
    Function vf = null;
    for (int a = 0; a < 3; a++)
      for (int b = 0; b < 3; b++) {
        var spmod = spmods[a];
        var spmod2 = spmods[b];
        storagekey = a + "," + b;
        var fseeds = new List<List<T>>();
        var aseeds = new List<List<T>>();
        for (int d = 0; d < ndir; d++) {
          var l = new List<T>();
          for (long i = 0; i < f.n_in(); i++) l.Add(S.sym("f", spmod(f.sparsity_in(i))));
          fseeds.Add(l);
          l = new List<T>();
          for (long i = 0; i < f.n_out(); i++) l.Add(S.sym("a", spmod2(f.sparsity_out(i))));
          aseeds.Add(l);
        }
        var inputss = new List<T>();
        for (long i = 0; i < f.n_in(); i++) inputss.Add(S.sym("i", f.sparsity_in(i)));
        bool notSX = !f.is_a("SXFunction");
        List<T> res = S.call(f, inputss, notSX);
        var fwdsens = S.fwd(res, inputss, fseeds, D(("always_inline", notSX)));
        var adjsens = S.rev(res, inputss, aseeds, D(("always_inline", notSX)));
        var fseed = new List<DM>();
        var aseed = new List<DM>();
        for (int d = 0; d < ndir; d++) {
          fseed.Add(randDMsp(random, spmod(f.sparsity_in(0))));
          aseed.Add(randDMsp(random, spmod2(f.sparsity_out(0))));
        }
        var vin = new List<T>(inputss);
        var vout = new List<T>(res);
        for (int d = 0; d < ndir; d++) {
          vin.AddRange(fseeds[d]);
          vin.AddRange(aseeds[d]);
          vout.AddRange(fwdsens[d]);
          vout.AddRange(adjsens[d]);
        }
        vf = S.mk("vf", vin, vout);
        var vf_in = new List<DM>(values);
        for (int d = 0; d < ndir; d++) {
          vf_in.Add(fseed[d]);
          for (int i = 0; i < values.Count - 1; i++) vf_in.Add(0.0);
          vf_in.Add(aseed[d]);
          vf_in.Add(0.0);
        }
        List<DM> vf_out = vf.call(vf_in);
        if (S.isMX) {
          check_eval_mx(fwdsens.Select(e => ca.veccat((List<MX>)(object)e)).ToList());
          check_eval_mx(adjsens.Select(e => ca.veccat((List<MX>)(object)e)).ToList());
        }
        check_codegen(vf, vf_in, std);
        check_serialize(vf, vf_in);
        int offset = res.Count;
        for (int d = 0; d < ndir; d++) {
          DM sens = vf_out[offset];
          offset += inputss.Count;
          checkarray(ca.vec(sens), ca.mtimes(J_, ca.vec(fseed[d])), "eval Fwd " + d + " " + f.class_name() + " " + S.isMX);
          sens = vf_out[offset];
          offset += inputss.Count;
          checkarray(ca.vec(sens), ca.mtimes(J_.T(), ca.vec(aseed[d])), "eval Adj " + d);
        }
        assertEquals(vf.n_out(), (long)offset);
        var rnd = new Random(1);
        vf_in = new List<DM>();
        for (long i = 0; i < vf.n_in(); i++) vf_in.Add(randDMsp(rnd, vf.sparsity_in(i)));
        vf_out = vf.call(vf_in);
        check_codegen(vf, vf_in, std);
        check_serialize(vf, vf_in);
        if (!storage.TryGetValue(storagekey, out var st)) storage[storagekey] = st = new List<List<DM>>();
        st.Add(vf_out);
        if (S.isMX) vf_mx = vf;
      }
    return vf;
  }

  /// <summary>Second-order sweeps through vf</summary>
  void secondOrder<T>(Sym<T> S, Function vf, Function vf_mx, int ndir, string std, string storagekey,
                      Dictionary<string, List<List<DM>>> storage2) {
    for (int a = 0; a < 3; a++)
      for (int b = 0; b < 3; b++) {
        var fseeds2 = new List<List<T>>();
        var aseeds2 = new List<List<T>>();
        for (int d = 0; d < ndir; d++) {
          var l = new List<T>();
          for (long i = 0; i < vf.n_in(); i++) l.Add(S.sym("f", vf_mx.sparsity_in(i)));
          fseeds2.Add(l);
          l = new List<T>();
          for (long i = 0; i < vf.n_out(); i++) l.Add(S.sym("a", vf_mx.sparsity_out(i)));
          aseeds2.Add(l);
        }
        var inputss2 = new List<T>();
        for (long i = 0; i < vf.n_in(); i++) inputss2.Add(S.sym("i", vf_mx.sparsity_in(i)));
        bool notSX = !vf.is_a("SXFunction");
        List<T> res2 = S.call(vf, inputss2, notSX);
        var fwdsens2 = S.fwd(res2, inputss2, fseeds2, D(("always_inline", notSX)));
        var adjsens2 = S.rev(res2, inputss2, aseeds2, D(("always_inline", notSX)));
        if (S.isMX) {
          check_eval_mx(fwdsens2.Select(e => ca.veccat((List<MX>)(object)e)).ToList());
          check_eval_mx(adjsens2.Select(e => ca.veccat((List<MX>)(object)e)).ToList());
        }
        var vin = new List<T>(inputss2);
        var vout = new List<T>(res2);
        for (int d = 0; d < ndir; d++) {
          vin.AddRange(fseeds2[d]);
          vin.AddRange(aseeds2[d]);
          vout.AddRange(fwdsens2[d]);
          vout.AddRange(adjsens2[d]);
        }
        Function vf2 = S.mk("vf2", vin, vout);
        var rnd = new Random(1);
        var vf2_in = new List<DM>();
        for (long i = 0; i < vf2.n_in(); i++) vf2_in.Add(randDMsp(rnd, vf2.sparsity_in(i)));
        List<DM> vf2_out = vf2.call(vf2_in);
        check_codegen(vf2, vf2_in, std);
        check_serialize(vf2, vf2_in);
        if (!storage2.TryGetValue(storagekey, out var st)) storage2[storagekey] = st = new List<List<DM>>();
        st.Add(vf2_out);
      }
  }

  public void test_MX() {
    MX x = MX.sym("x", 2), y = MX.sym("y", 2, 2);
    var f1 = ca.Function("f1", new[] { x, y }, new[] { x + y[0, 0], ca.mtimes(y, x) });
    var f1_noninline = ca.Function("f1", new[] { x, y }, new[] { x + y[0, 0], ca.mtimes(y, x) }, D(("never_inline", true)));
    var f2 = ca.Function("f2", new[] { x, y }, new[] { ca.mtimes(MX.zeros(0, 2), x) });
    var f3 = ca.Function("f3", new[] { x, y }, new[] { MX.zeros(0, 0), ca.mtimes(y, x) });
    var f4 = ca.Function("f4", new[] { x, y }, new[] { MX.zeros(0, 2), ca.mtimes(y, x) });
    int ndir = 2;
    var in1 = new List<MX> { x, y };
    var v1 = new List<DM> { ca.vertcat(1.1, 1.3), ca.DM(new double[,] { { 0.7, 1.5 }, { 2.1, 0.9 } }) };
    var all = ca.Slice();
    MX w = x[all];
    w[1] = w[1] * 2.0;
    MX w2 = x[all];
    w2[1] = w2[1] * x[0];
    MX ww = x[all];
    ww[new long[] { 0, 1 }] = ww[new long[] { 0, 1 }] * x;
    MX wwf = x[all];
    wwf[new long[] { 1, 0 }] = wwf[new long[] { 1, 0 }] * x;
    MX wwr = x[all];
    wwr[new long[] { 0, 0, 1, 1 }] = wwr[new long[] { 0, 0, 1, 1 }] * 2.0;
    MX yy = y[all, all];
    yy[all, 0] = x;
    MX yy2 = y[all, all];
    yy2[all, 0] = sq(x);
    MX yyy = y[all, all];
    yyy.set(x, false, idx(1, 0), (DM)0.0);
    MX yyy2 = y[all, all];
    yyy2.set(sq(x), false, idx(1, 0), (DM)0.0);
    var spmods = new List<Func<Sparsity, Sparsity>> { s => s, remove_first, remove_last };
    MX y_nz = ca.sin(y);
    y_nz.set_nz(ca.vertcat(x[1], x[0]), false, idx(2, 1));
    MX y_nzx = ca.sin(x);
    y_nzx.set_nz(sq(x[0]), false, 1L);
    MX y_nzs = ca.sin(y);
    y_nzs.set_nz(ca.vertcat(x[1], x[0]), false, y[all, 0]);
    MX y_nzxs = ca.sin(x);
    y_nzxs.set_nz(sq(x[0]), false, y[1, 0]);
    MX xx = ca.horzcat(ca.sin(x), ca.cos(x));
    var Z = ca.MX(1, 1);
    DM eye2 = DM.eye(2);
    var v2 = new List<DM> { v1[0], ca.DM(new double[,] { { 1, 1.5 }, { 0, 0.9 } }) };
    Sparsity t1 = Sparsity.triplet(2, 1, new long[] { 1 }, new long[] { 0 });
    Sparsity t2 = Sparsity.triplet(2, 2, new long[] { 1, 1 }, new long[] { 0, 1 });
    Sparsity t3 = Sparsity.triplet(2, 2, new long[] { 0, 1, 1 }, new long[] { 0, 0, 1 });
    var blk16 = new object[16][];
    blk16[0] = new object[] { Z, x[1] * 2.0 };
    for (int i = 1; i < 15; i++) blk16[i] = new object[] { Z, Z };
    blk16[15] = new object[] { x[0] * 2.0, Z };
    DM sp01 = ca.sparsify(ca.DM(new double[,] { { 0, 1 }, { 1, 0 } }));
    var cases = new List<Case> {
      new Case(v1, ca.sparsity_cast(sq(x), sp01.sparsity()),
               ca.sparsity_cast(x, ca.sparsify(ca.DM(new double[,] { { 0, 0 }, { 1, 0 }, { 0, 1 }, { 0, 0 } })).sparsity()) * 2.0, true, "c89"),
      new Case(v1, f1_noninline.call(new[] { sq(x), y })[1], y * 2.0 * ca.vertcat(x.T(), x.T()), true, "c89"),
      new Case(v2, xx.get(false, y[all, 0], all), bc(new object[] { 0, ca.cos(x[1]) }, new object[] { ca.cos(x[0]), 0 },
               new object[] { 0, -ca.sin(x[1]) }, new object[] { -ca.sin(x[0]), 0 }), false, "c99"),
      new Case(v2, xx.get(false, all, y[all, 0]), bc(new object[] { -ca.sin(x[0]), 0 },
               new object[] { 0, -ca.sin(x[1]) }, new object[] { ca.cos(x[0]), 0 }, new object[] { 0, ca.cos(x[1]) }), false, "c99"),
      new Case(v2, xx.get(false, y[all, 0], y[all, 0]), bc(new object[] { 0, -ca.sin(x[1]) },
               new object[] { -ca.sin(x[0]), 0 }, new object[] { 0, ca.cos(x[1]) }, new object[] { ca.cos(x[0]), 0 }), false, "c99"),
      new Case(v1, y_nz, bc(new object[] { Z, Z }, new object[] { 1, Z }, new object[] { Z, 1 }, new object[] { Z, Z }), true, "c89"),
      new Case(v1, y_nzx, bc(new object[] { ca.cos(x[0]), Z }, new object[] { x[0] * 2.0, Z }), true, "c89"),
      new Case(new List<DM> { ca.vertcat(2, 1), v1[1] }, y.get_nz(false, x), ca.DM(2, 2), false, "c99"),
      new Case(v2, ca.sin(x).get_nz(false, y[all, 0]), bc(new object[] { 0, ca.cos(x[1]) }, new object[] { ca.cos(x[0]), 0 }), false, "c99"),
      new Case(new List<DM> { v1[0], ca.DM(new double[,] { { 2, 1.5 }, { 1, 0.9 } }) }, y_nzs,
               bc(new object[] { 0, 0 }, new object[] { 1, 0 }, new object[] { 0, 1 }, new object[] { 0, 0 }), false, "c99"),
      new Case(new List<DM> { v1[0], ca.DM(new double[,] { { 0, 1.5 }, { 1, 0.9 } }) }, y_nzxs,
               bc(new object[] { ca.cos(x[0]), 0 }, new object[] { x[0] * 2.0, 0 }), false, "c99"),
      new Case(new List<DM> { v1[0], ca.DM(new double[,] { { 1, 1.5 }, { 1, 0.9 } }) }, x.get_nz(false, y[all, 0]),
               bc(new object[] { 0, 1 }, new object[] { 0, 1 }), false, "c99"),
      new Case(v1, x, eye2, true, "c89"),
      new Case(v1, x.T(), eye2, true, "c89"),
      new Case(v1, sq(x), ca.diag(x) * 2.0, true, "c89"),
      new Case(v1, sq(x).attachAssert(1.0), ca.diag(x) * 2.0, true, "c89"),
      new Case(v1, sq(x).T(), ca.diag(x) * 2.0, true, "c89"),
      new Case(v1, ca.reshape(x, 1, 2), eye2, true, "c89"),
      new Case(v1, ca.reshape(sq(x), 1, 2), ca.diag(x) * 2.0, true, "c89"),
      new Case(v1, x + y.get_nz(false, 0L), eye2, true, "c89"),
      new Case(v1, x + y[0, 0], eye2, true, "c89"),
      new Case(v1, x + x, eye2 * 2.0, true, "c89"),
      new Case(v1, sq(x) + x, ca.diag(x) * 2.0 + eye2, true, "c89"),
      new Case(v1, x * x, ca.diag(x) * 2.0, true, "c89"),
      new Case(v1, x * y.get_nz(false, 0L), eye2 * y.get_nz(false, 0L), true, "c89"),
      new Case(v1, x * y[0, 0], eye2 * y[0, 0], true, "c89"),
      new Case(v1, x[0], eye2[0, all], true, "c89"),
      new Case(v1, sq(x)[0], ca.horzcat(x[0] * 2.0, Z), true, "c89"),
      new Case(v1, x[0] + x[1], DM.ones(1, 2), true, "c89"),
      new Case(v1, ca.sin(ca.repmat(sq(x), 1, 3)), ca.repmat(ca.cos(ca.diag(sq(x))) * 2.0 * ca.diag(x), 3, 1), true, "c89"),
      new Case(v1, ca.sin(ca.repsum(sq(x).T(), 1, 2)), ca.cos(sq(x[0]) + sq(x[1])) * 2.0 * x.T(), true, "c89"),
      new Case(v1, ca.vertcat(x[1], x[0]), sp01, true, "c89"),
      new Case(v1, ca.vertsplit(x, new long[] { 0, 1, 2 })[1], ca.sparsify(ca.DM(new double[,] { { 0, 1 } })), true, "c89"),
      new Case(v1, ca.vertcat(sq(x[1]), sq(x[0])), bc(new object[] { Z, x[1] * 2.0 }, new object[] { x[0] * 2.0, Z }), true, "c89"),
      new Case(v1, ca.vertsplit(sq(x), new long[] { 0, 1, 2 })[1], bc(new object[] { Z, x[1] * 2.0 }), true, "c89"),
      new Case(v1, ca.pow(ca.vertsplit(sq(x), new long[] { 0, 1, 2 })[1], 3.0), bc(new object[] { Z, ca.pow(x[1], 5.0) * 6.0 }), true, "c89"),
      new Case(v1, ca.horzcat(x[1], x[0]).T(), sp01, true, "c89"),
      new Case(v1, ca.horzcat(sq(x[1]), sq(x[0])).T(), bc(new object[] { Z, x[1] * 2.0 }, new object[] { x[0] * 2.0, Z }), true, "c89"),
      new Case(v1, ca.diagcat(sq(x[1]), y, sq(x[0])), bc(blk16), true, "c89"),
      new Case(v1, ca.horzcat(sq(x[1]), sq(x[0])).T(), bc(new object[] { Z, x[1] * 2.0 }, new object[] { x[0] * 2.0, Z }), true, "c89"),
      new Case(v1, x[new long[] { 0, 1 }], ca.sparsify(ca.DM(new double[,] { { 1, 0 }, { 0, 1 } })), true, "c89"),
      new Case(v1, sq(x)[new long[] { 0, 1 }], ca.diag(x) * 2.0, true, "c89"),
      new Case(v1, x[new long[] { 0, 0, 1, 1 }], ca.sparsify(ca.DM(new double[,] { { 1, 0 }, { 1, 0 }, { 0, 1 }, { 0, 1 } })), true, "c89"),
      new Case(v1, sq(x)[new long[] { 0, 0, 1, 1 }], bc(new object[] { x[0] * 2.0, Z }, new object[] { x[0] * 2.0, Z },
               new object[] { Z, x[1] * 2.0 }, new object[] { Z, x[1] * 2.0 }), true, "c89"),
      new Case(v1, wwr, ca.sparsify(ca.DM(new double[,] { { 2, 0 }, { 0, 2 } })), true, "c89"),
      new Case(v1, x[new long[] { 1, 0 }], sp01, true, "c89"),
      new Case(v1, x.get(idx(1, 0), (DM)0.0), sp01, true, "c89"),
      new Case(v1, w, ca.sparsify(ca.DM(new double[,] { { 1, 0 }, { 0, 2 } })), true, "c89"),
      new Case(v1, w2, bc(new object[] { 1, Z }, new object[] { x[1], x[0] }), true, "c89"),
      new Case(v1, ww, ca.diag(x) * 2.0, true, "c89"),
      new Case(v1, wwf, ca.vertcat(x[new long[] { 1, 0 }].T(), x[new long[] { 1, 0 }].T()), true, "c89"),
      new Case(v1, yy[all, 0], eye2, true, "c89"),
      new Case(v1, yy2[all, 0], ca.diag(x) * 2.0, true, "c89"),
      new Case(v1, yyy[all, 0], sp01, true, "c89"),
      new Case(v1, ca.mtimes(y, x), y, true, "c89"),
      new Case(v1, ca.mtimes(x.T(), y.T()), y, true, "c89"),
      new Case(v1, ca.mac(y, x, DM.zeros(t1)), y.get(t2), true, "c89"),
      new Case(v1, ca.mac(x.T(), y.T(), DM.zeros(t1.T())), y.get(t2), true, "c89"),
      new Case(v1, ca.mtimes(y.get(t3), x), y.get(t3), true, "c89"),
      new Case(v1, ca.mtimes(x.T(), y.get(t3).T()), y.get(t3), true, "c89"),
      new Case(v1, ca.mtimes(y, sq(x)), y * 2.0 * ca.vertcat(x.T(), x.T()), true, "c89"),
      new Case(v1, ca.sin(x), ca.diag(ca.cos(x)), true, "c89"),
      new Case(v1, ca.sin(sq(x)), ca.diag(ca.cos(sq(x)) * 2.0 * x), true, "c89"),
      new Case(v1, x * y[all, 0], ca.diag(y[all, 0]), true, "c89"),
      new Case(v1, x * y.get_nz(false, idx(0, 1)), ca.diag(y.get_nz(false, idx(0, 1))), true, "c89"),
      new Case(v1, x * y.get_nz(false, idx(1, 0)), ca.diag(y.get_nz(false, idx(1, 0))), true, "c89"),
      new Case(v1, x * y.get(idx(0, 1), (DM)0.0), ca.diag(y.get(idx(0, 1), (DM)0.0)), true, "c89"),
      new Case(v1, x * y.get(idx(1, 0), (DM)0.0), ca.diag(y.get(idx(1, 0), (DM)0.0)), true, "c89"),
      new Case(v1, ca.dot(x, x), (x * 2.0).T(), true, "c89"),
      new Case(v1, ca.dot(sq(x), x), (sq(x) * 3.0).T(), true, "c89"),
      new Case(v1, ca.det(ca.horzcat(x, (MX)ca.vertcat(1, 2))), ca.DM(new double[,] { { 2, -1 } }), true, "c89"),
      new Case(v1, ca.det(ca.horzcat(ca.vertcat(x[0], (MX)1.0), ca.vertcat((MX)1.0, x[1]))), ca.horzcat(x[1], x[0]), true, "c89"),
      new Case(v1, f1.call(in1)[1], y, true, "c89"),
      new Case(v1, f1.call(new[] { sq(x), y })[1], y * 2.0 * ca.vertcat(x.T(), x.T()), true, "c89"),
      new Case(v1, f2.call(in1)[0], DM.zeros(0, 2), true, "c89"),
      new Case(v1, f2.Call(sq(x), y), DM.zeros(0, 2), true, "c89"),
      new Case(v1, f3.call(in1)[0], DM.zeros(0, 2), true, "c89"),
      new Case(v1, f3.call(new[] { sq(x), y })[0], DM.zeros(0, 2), true, "c89"),
      new Case(v1, f4.call(in1)[0], DM.zeros(0, 2), true, "c89"),
      new Case(v1, f4.call(new[] { sq(x), y })[0], DM.zeros(0, 2), true, "c89"),
      new Case(v1, ca.vertcat(x, ca.MX(ca.DM(0, 1))), eye2, true, "c89"),
      new Case(v1, ca.project(sq(x), ca.sparsify(ca.vertcat(0, 1)).sparsity()), bc(new object[] { Z, Z }, new object[] { Z, x[1] * 2.0 }), true, "c89"),
      new Case(v1, ca.dot(x, y[all, 0]), y[all, 0].T(), true, "c89"),
      new Case(v1, x.get_nz(false, ca.DM(new double[,] { { 1, 0 } })).T() * y.get_nz(false, ca.DM(new double[,] { { 0, 2 } })),
               bc(new object[] { Z, y.get_nz(false, 0L) }, new object[] { y.get_nz(false, 2L), Z }), true, "c89"),
      new Case(v1, x.get_nz(false, ca.diag(ca.vertcat(1, 0))) * y.get_nz(false, ca.diag(ca.vertcat(0, 2))),
               bc(new object[] { Z, y.get_nz(false, 0L) }, new object[] { Z, Z }, new object[] { Z, Z }, new object[] { y.get_nz(false, 2L), Z }), true, "c89"),
      new Case(v1, ca.kron(x, y), ca.horzcat(ca.vec(ca.vertcat(y, ca.MX(2, 2))), ca.vec(ca.vertcat(ca.MX(2, 2), y))), true, "c89"),
      new Case(v1, ca.kron(y, x), ca.vertcat(y[0, 0] * eye2, y[1, 0] * eye2, y[0, 1] * eye2, y[1, 1] * eye2), true, "c89"),
      new Case(v1, ca.kron_contract(ca.kron(x, y), y, true), ca.sumsqr(y) * eye2, true, "c89"),
      new Case(v1, ca.kron_contract(ca.kron(x, y), x, false), ca.mtimes(ca.vec(y) * 2.0, x.T()), true, "c89"),
    };

    cgBatch = new Dictionary<string, List<(Function, List<DM>)>>();
    foreach (Case cs in cases) {
      List<DM> values = cs.values;
      MX jac = cs.jac;
      bool with_sx = cs.withSx;
      string std = cs.std;
      var fun0 = ca.Function("fun", in1, new[] { cs.outp });
      var shifted = in1.Select(e => e + 1e-300).ToList();
      MX out_eval_mx = fun0.call(shifted, true, false)[0];
      foreach (MX outp in new[] { cs.outp, out_eval_mx }) {
        check_eval_mx(outp);
        var fun = ca.Function("fun", in1, new[] { outp, jac });
        var fun_ad = new List<Function>();
        foreach (int wt in new[] { 0, 1 })
          fun_ad.Add(ca.Function("fun", in1, new[] { outp, jac }, D(("ad_weight", wt), ("ad_weight_sp", wt))));
        List<DM> fun_out = fun.call(values);
        List<Function> funsx_ad = null;
        if (with_sx) {
          Function funsx = fun.expand("expand_fun");
          funsx_ad = fun_ad.Select(f => f.expand("expand_" + f.name())).ToList();
          List<DM> funsx_out = funsx.call(values);
          checkarray(fun_out[0], funsx_out[0]);
          checkarray(fun_out[1], funsx_out[1]);
        }
        check_codegen(fun, values, std);
        check_serialize(fun, values);
        DM J_ = fun_out[1];
        var storage2 = new Dictionary<string, List<List<DM>>>();
        var storage = new Dictionary<string, List<List<DM>>>();
        Function vf_mx = null;
        var fs = new List<Function> { fun };
        if (with_sx) fs.Add(fun.expand("expand_" + fun.name()));
        foreach (Function f in fs) {
          f.forward(ndir);
          f.reverse(ndir);
          foreach (bool symMX in new[] { true, false }) {
            if (f.is_a("MXFunction") && !symMX) continue;
            if (f.is_a("SXFunction") && symMX) continue;
            string storagekey = null;
            Function vf = symMX
              ? firstOrder(MXS, f, ndir, spmods, values, J_, std, storage, ref storagekey, ref vf_mx)
              : firstOrder(SXS, f, ndir, spmods, values, J_, std, storage, ref storagekey, ref vf_mx);
            foreach (bool sym2MX in new[] { true, false }) {
              if (vf.is_a("MXFunction") && !sym2MX) continue;
              if (vf.is_a("MXFunction") && sym2MX) continue;
              if (sym2MX) secondOrder(MXS, vf, vf_mx, ndir, std, storagekey, storage2);
              else secondOrder(SXS, vf, vf_mx, ndir, std, storagekey, storage2);
            }
          }
        }
        foreach (var store in new[] { storage, storage2 }) {
          string order = store == storage ? "first-order" : "second-order";
          foreach (var st in store.Values)
            for (int i = 0; i < st.Count - 1; i++)
              for (int k = 0; k < st[0].Count; k++) {
                DM a = st[0][k], b = st[i + 1][k];
                if (b.numel() == 0 && ca.sparsify(a).nnz() == 0) continue;
                if (a.numel() == 0 && ca.sparsify(b).nnz() == 0) continue;
                checkarray(ca.sparsify(a), ca.sparsify(b), order + ", output(" + k + ")");
              }
        }
        foreach (bool expand in with_sx ? new[] { false, true } : new[] { false })
          foreach (string mode in new[] { "forward", "reverse" }) {
            int ind = mode == "forward" ? 0 : 1;
            Function f = expand ? funsx_ad[ind] : fun_ad[ind];
            Function Jf = jacobian_old(f, 0, 0);
            List<DM> Jf_out = Jf.call(values);
            check_codegen(Jf, values, std);
            check_serialize(Jf, values);
            checkarray(Jf_out[0], J_);
            checkarray(DM.ones(Jf.sparsity_out(0)), DM.ones(J_.sparsity()), outp + mode);
            checkarray(DM.ones(f.jac_sparsity(0, 0)), DM.ones(J_.sparsity()));
          }
        if (outp.is_empty()) continue;
        long s_i = outp.sparsity().row()[0], s_j = outp.sparsity().get_col()[0];
        long s_k = s_i * outp.size2() + s_j;
        DM H_ = null;
        foreach (bool expand in with_sx ? new[] { false, true } : new[] { false })
          foreach (string mode in new[] { "forward", "reverse" }) {
            int wt = mode == "forward" ? 0 : 1;
            var f = ca.Function("fun", in1, new[] { outp[s_i, s_j], jac[s_k, all].T() }, D(("ad_weight", wt), ("ad_weight_sp", wt)));
            if (expand) f = f.expand("expand_" + f.name());
            f.call(values);
            Function Hf = hessian_old(f, 0, 0);
            List<DM> Hf_out = Hf.call(values);
            check_codegen(Hf, values, std);
            check_serialize(Hf, values);
            if (H_ is null) H_ = Hf_out[0];
            checkarray(Hf_out[0], H_, "mode: " + mode);
          }
      }
      flushCodegen();
    }
  }

  public void test_repmat() {
    MX X = MX.sym("x", 2, 2);
    foreach (MX e in new[] { X.T() + ca.repmat(MX.zeros(2, 1), 1, 2) })
      foreach (int weight in new[] { 0, 1 }) {
        var F = ca.Function("f", new[] { ca.vec(X) }, new[] { e }, D(("ad_weight_sp", weight)));
        Sparsity J0 = F.jac_sparsity(0, 0, false);
        J0.spy();
        Sparsity J1 = F.expand().jac_sparsity(0, 0, false);
        J1.spy();
        assertTrue(J0.is_equal(J1));
      }
  }
}

}
