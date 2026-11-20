package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.UnaryOperator;

import org.casadi.CodeGenerator;
import org.casadi.DM;
import org.casadi.GenericMatrix;
import org.casadi.GlobalOptions;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.Sparsity;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

/** Port of test/python/ad.py */
public class ADTest extends CasadiTestCase {

  private static final String[] SHAPES = {"column", "row", "matrix"};
  private static final String[] TYPES = {"dense", "sparse"};
  private static final double CODEGEN_TOL = 0.5e-15;

  private Map<String, SX> sxinputs, sxoutputs;
  private Map<String, MX> mxinputs;
  private Map<String, UnaryOperator<MX>> mxoutputs;

  private static MX nz(MX x, long k) {
    return x.get_nz(false, k);
  }

  private static MX temp1(MX xyz) {
    MX X = new MX(6, 1);
    X.set(0, 0, nz(xyz, 0));
    X.set(2, 0, nz(xyz, 0).plus(ca.times(2.0, nz(xyz, 1).power(2.0))));
    X.set(4, 0, nz(xyz, 0).plus(ca.times(2.0, nz(xyz, 1).power(3.0))).plus(ca.times(3.0, nz(xyz, 2).power(4.0))));
    X.set(5, 0, nz(xyz, 3));
    return X;
  }

  private static MX temp2(MX xyz) {
    MX X = new MX(1, 6);
    X.set(0, 0, nz(xyz, 0));
    X.set(0, 2, nz(xyz, 0).plus(ca.times(2.0, nz(xyz, 1).power(2.0))));
    X.set(0, 4, nz(xyz, 0).plus(ca.times(2.0, nz(xyz, 1).power(3.0))).plus(ca.times(3.0, nz(xyz, 2).power(4.0))));
    X.set(0, 5, nz(xyz, 3));
    return X;
  }

  private static List<MX> outList(MX xyz) {
    return List.of(nz(xyz, 0), nz(xyz, 0).plus(ca.times(2.0, nz(xyz, 1).power(2.0))),
        nz(xyz, 0).plus(ca.times(2.0, nz(xyz, 1).power(3.0))).plus(ca.times(3.0, nz(xyz, 2).power(4.0))), nz(xyz, 3));
  }

  @BeforeEach
  public void setUp() {
    SX x = SX.sym("x"), y = SX.sym("y"), z = SX.sym("z"), w = SX.sym("w");
    SX out = new SX(6, 1);
    out.set(0, 0, x);
    out.set(2, 0, x.plus(ca.times(2.0, y.power(2.0))));
    out.set(4, 0, x.plus(ca.times(2.0, y.power(3.0))).plus(ca.times(3.0, z.power(4.0))));
    out.set(5, 0, w);
    SX inp = new SX(6, 1);
    inp.set(0, 0, x);
    inp.set(2, 0, y);
    inp.set(4, 0, z);
    inp.set(5, 0, w);
    Sparsity sp = new Sparsity(1, 6, new long[] {0, 1, 1, 2, 2, 3, 4}, new long[] {0, 0, 0, 0}).T();
    Sparsity spT = new Sparsity(6, 1, new long[] {0, 4}, new long[] {0, 2, 4, 5}).T();
    SX xyzw = ca.vertcat(List.of(x, y, z, w));
    sxinputs = new LinkedHashMap<>();
    sxinputs.put("column dense", xyzw);
    sxinputs.put("column sparse", inp);
    sxinputs.put("row dense", xyzw.T());
    sxinputs.put("row sparse", inp.T());
    sxinputs.put("matrix dense", ca.reshape(xyzw, 2, 2));
    sxinputs.put("matrix sparse", ca.reshape(inp, 3, 2));
    mxinputs = new LinkedHashMap<>();
    mxinputs.put("column dense", MX.sym("xyzw", 4, 1));
    mxinputs.put("column sparse", MX.sym("xyzw", sp));
    mxinputs.put("row dense", MX.sym("xyzw", 1, 4));
    mxinputs.put("row sparse", MX.sym("xyzw", spT));
    mxinputs.put("matrix dense", MX.sym("xyzw", 2, 2));
    mxinputs.put("matrix sparse", MX.sym("xyzw", ca.reshape(inp, 3, 2).sparsity()));
    mxoutputs = new LinkedHashMap<>();
    mxoutputs.put("column dense", xyz -> ca.vertcat(outList(xyz)));
    mxoutputs.put("column sparse", ADTest::temp1);
    mxoutputs.put("row dense", xyz -> ca.horzcat(outList(xyz)));
    mxoutputs.put("row sparse", ADTest::temp2);
    mxoutputs.put("matrix dense", xyz -> ca.reshape(ca.vertcat(outList(xyz)), 2, 2));
    mxoutputs.put("matrix sparse", xyz -> ca.reshape(temp1(xyz), 3, 2));
    SX o = ca.vertcat(List.of(x, x.plus(ca.times(2.0, y.power(2.0))),
        x.plus(ca.times(2.0, y.power(3.0))).plus(ca.times(3.0, z.power(4.0))), w));
    sxoutputs = new LinkedHashMap<>();
    sxoutputs.put("column dense", o);
    sxoutputs.put("column sparse", out);
    sxoutputs.put("row dense", o.T());
    sxoutputs.put("row sparse", out.T());
    sxoutputs.put("matrix dense", ca.reshape(o, 2, 2));
    sxoutputs.put("matrix sparse", ca.reshape(out, 3, 2));
  }

  private static DM jacobians(String inputtype, String outputtype, double[] n) {
    double y = n[1], z = n[2];
    String k = inputtype + " " + outputtype;
    double[][] r;
    if (k.equals("dense dense"))
      r = new double[][] {{1, 0, 0, 0}, {1, 4 * y, 0, 0}, {1, 6 * y * y, 12 * z * z * z, 0}, {0, 0, 0, 1}};
    else if (k.equals("dense sparse"))
      r = new double[][] {{1, 0, 0, 0}, {0, 0, 0, 0}, {1, 4 * y, 0, 0}, {0, 0, 0, 0}, {1, 6 * y * y, 12 * z * z * z, 0}, {0, 0, 0, 1}};
    else if (k.equals("sparse dense"))
      r = new double[][] {{1, 0, 0, 0, 0, 0}, {1, 0, 4 * y, 0, 0, 0}, {1, 0, 6 * y * y, 0, 12 * z * z * z, 0}, {0, 0, 0, 0, 0, 1}};
    else
      r = new double[][] {{1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {1, 0, 4 * y, 0, 0, 0}, {0, 0, 0, 0, 0, 0},
          {1, 0, 6 * y * y, 0, 12 * z * z * z, 0}, {0, 0, 0, 0, 0, 1}};
    return new DM(r);
  }

  private static DM nonzeroIndicator(DM J) {
    double[][] a = J.full();
    for (double[] row : a)
      for (int j = 0; j < row.length; j++) row[j] = row[j] != 0 ? 1 : 0;
    return new DM(a);
  }

  private static org.casadi.Function jacobian_old(org.casadi.Function f, long i, long j) {
    List<String> out = new ArrayList<>(List.of("jac:" + f.name_out(j) + ":" + f.name_in(i)));
    out.addAll(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), out);
  }

  private static org.casadi.Function hessian_old(org.casadi.Function f, long i, long j) {
    List<String> out = new ArrayList<>(List.of("hess:" + f.name_out(j) + ":" + f.name_in(i) + ":" + f.name_in(i),
        "grad:" + f.name_out(j) + ":" + f.name_in(i)));
    out.addAll(f.name_out());
    return f.factory(f.name() + "_hess", f.name_in(), out);
  }

  private static List<List<DM>> seedLists(Sparsity sp, double[][] seeds) {
    List<List<DM>> r = new ArrayList<>();
    for (double[] s : seeds) r.add(List.of(new DM(sp, new DM(s))));
    return r;
  }

  private static <T> List<List<T>> wrapSeeds(List<List<DM>> s, Function<DM, T> conv) {
    List<List<T>> r = new ArrayList<>();
    for (List<DM> l : s) r.add(List.of(conv.apply(l.get(0))));
    return r;
  }

  private static final double[][] SEEDS = {{1, 0, 0, 0}, {0, 2, 0, 0}, {1.2, 4.8, 7.9, 4.6}};

  @Test
  public void test_SXeval_sx() {
    double[] n = {1.2, 2.3, 7, 1.4};
    for (String inputshape : SHAPES)
      for (String outputshape : SHAPES)
        for (String inputtype : TYPES)
          for (String outputtype : TYPES) {
            message("eval_sx on SX. Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
            org.casadi.Function f = new org.casadi.Function("f", List.of(sxinputs.get(inputshape + " " + inputtype)),
                List.of(sxoutputs.get(outputshape + " " + outputtype)));
            DM f_in = new DM(f.sparsity_in(0), new DM(n));
            DM r = (DM) f.invoke(f_in);
            DM J = jacobians(inputtype, outputtype, n);
            SX y = SX.sym("y", f.sparsity_in(0));
            List<List<DM>> fseeds = seedLists(f.sparsity_in(0), SEEDS), aseeds = seedLists(f.sparsity_out(0), SEEDS);
            SX res = (SX) f.invoke(y);
            List<List<SX>> fwdsens = ca.forward(List.of(res), List.of(y), wrapSeeds(fseeds, SX::new));
            List<List<SX>> adjsens = ca.reverse(List.of(res), List.of(y), wrapSeeds(aseeds, SX::new));
            org.casadi.Function fe = new org.casadi.Function("fe", List.of(y), List.of(res));
            DM re = (DM) fe.invoke(f_in);
            checkarray(r, re);
            for (int i = 0; i < fwdsens.size(); i++) {
              fe = new org.casadi.Function("fe", List.of(y), List.of(fwdsens.get(i).get(0)));
              re = (DM) fe.invoke(f_in);
              checkarray(ca.vec(re), ca.mtimes(J, ca.vec(fseeds.get(i).get(0))), "AD");
            }
            for (int i = 0; i < adjsens.size(); i++) {
              fe = new org.casadi.Function("fe", List.of(y), List.of(adjsens.get(i).get(0)));
              re = (DM) fe.invoke(f_in);
              checkarray(ca.vec(re), ca.mtimes(J.T(), ca.vec(aseeds.get(i).get(0))), "AD");
            }
          }
  }

  @Test
  public void test_MXeval_mx() {
    double[] n = {1.2, 2.3, 7, 1.4};
    for (String inputshape : SHAPES)
      for (String outputshape : SHAPES)
        for (String inputtype : TYPES)
          for (String outputtype : TYPES) {
            message("eval_mx on MX. Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
            MX in = mxinputs.get(inputshape + " " + inputtype);
            org.casadi.Function f = new org.casadi.Function("f", List.of(in),
                List.of(mxoutputs.get(outputshape + " " + outputtype).apply(in)));
            DM f_in = new DM(f.sparsity_in(0), new DM(n));
            DM r = (DM) f.invoke(f_in);
            DM J = jacobians(inputtype, outputtype, n);
            MX y = MX.sym("y", f.sparsity_in(0));
            List<List<DM>> fseeds = seedLists(f.sparsity_in(0), SEEDS), aseeds = seedLists(f.sparsity_out(0), SEEDS);
            MX res = (MX) f.invoke(y);
            List<List<MX>> fwdsens = ca.forward(List.of(res), List.of(y), wrapSeeds(fseeds, MX::new));
            List<List<MX>> adjsens = ca.reverse(List.of(res), List.of(y), wrapSeeds(aseeds, MX::new));
            org.casadi.Function fe = new org.casadi.Function("fe", List.of(y), List.of(res));
            DM re = (DM) fe.invoke(f_in);
            checkarray(r, re);
            for (int i = 0; i < fwdsens.size(); i++) {
              fe = new org.casadi.Function("fe", List.of(y), List.of(fwdsens.get(i).get(0)));
              re = (DM) fe.invoke(f_in);
              checkarray(ca.vec(re), ca.mtimes(J, ca.vec(fseeds.get(i).get(0))), "AD");
            }
            for (int i = 0; i < adjsens.size(); i++) {
              fe = new org.casadi.Function("fe", List.of(y), List.of(adjsens.get(i).get(0)));
              re = (DM) fe.invoke(f_in);
              checkarray(ca.vec(re), ca.mtimes(J.T(), ca.vec(aseeds.get(i).get(0))), "AD");
            }
          }
  }

  @Test
  public void test_MXeval_sx() {
    Assumptions.assumeTrue(System.getenv("CASADI_JAVA_KNOWN_BUGS") != null, "known bug (not implemented)");
    double[] n = {1.2, 2.3, 7, 1.4};
    for (String inputshape : SHAPES)
      for (String outputshape : SHAPES)
        for (String inputtype : TYPES)
          for (String outputtype : TYPES) {
            message("eval_sx on MX. Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
            MX in = mxinputs.get(inputshape + " " + inputtype);
            org.casadi.Function f = new org.casadi.Function("f", List.of(in),
                List.of(mxoutputs.get(outputshape + " " + outputtype).apply(in)));
            List<DM> f_in = new ArrayList<>();
            for (long i = 0; i < f.n_in(); i++) f_in.add(new DM(0.0));
            f_in.set(0, new DM(n));
            DM r = f.call(f_in).get(0);
            DM J = jacobians(inputtype, outputtype, n);
            SX y = SX.sym("y", f.sparsity_in(0));
            List<List<DM>> fseeds = seedLists(f.sparsity_in(0), SEEDS), aseeds = seedLists(f.sparsity_out(0), SEEDS);
            SX res = (SX) f.invoke(y);
            List<List<SX>> fwdsens = ca.forward(List.of(res), List.of(y), wrapSeeds(fseeds, SX::new));
            List<List<SX>> adjsens = ca.reverse(List.of(res), List.of(y), wrapSeeds(aseeds, SX::new));
            org.casadi.Function fe = new org.casadi.Function("fe", List.of(y), List.of(res));
            List<DM> fe_in = List.of(new DM(n));
            checkarray(r, fe.call(fe_in).get(0));
            for (int i = 0; i < fwdsens.size(); i++) {
              fe = new org.casadi.Function("fe", List.of(y), List.of(fwdsens.get(i).get(0)));
              checkarray(ca.vec(fe.call(fe_in).get(0).T()), ca.mtimes(J, ca.vec(fseeds.get(i).get(0).T())), "AD");
            }
            for (int i = 0; i < adjsens.size(); i++) {
              fe = new org.casadi.Function("fe", List.of(y), List.of(adjsens.get(i).get(0)));
              checkarray(ca.vec(fe.call(fe_in).get(0).T()), ca.mtimes(J.T(), ca.vec(aseeds.get(i).get(0).T())), "AD");
            }
          }
  }

  @Test
  public void test_MXeval_sx_reduced() {
    double[] n = {1.2, 2.3, 7, 1.4};
    for (String inputshape : SHAPES)
      for (String outputshape : SHAPES)
        for (String inputtype : TYPES)
          for (String outputtype : TYPES) {
            message("eval_sx on MX. Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
            MX in = mxinputs.get(inputshape + " " + inputtype);
            org.casadi.Function f = new org.casadi.Function("f", List.of(in),
                List.of(mxoutputs.get(outputshape + " " + outputtype).apply(in)));
            DM f_in = new DM(f.sparsity_in(0), new DM(n));
            DM r = (DM) f.invoke(f_in);
            SX y = SX.sym("y", f.sparsity_in(0));
            SX res = (SX) f.invoke(y);
            ca.forward(List.of(res), List.of(y), new ArrayList<List<SX>>());
            ca.reverse(List.of(res), List.of(y), new ArrayList<List<SX>>());
            new org.casadi.Function("fe", List.of(y), List.of(res));
            DM re = (DM) f.invoke(f_in);
            checkarray(r, re);
          }
  }

  private static Map<String, Object> adOpts(String mode) {
    int w = mode.equals("forward") ? 0 : 1;
    return Map.of("ad_weight", w, "ad_weight_sp", w);
  }

  @Test
  public void test_Jacobian() {
    double[] n = {1.2, 2.3, 7, 4.6};
    for (String inputshape : SHAPES)
      for (String outputshape : SHAPES)
        for (String inputtype : TYPES)
          for (String outputtype : TYPES)
            for (String mode : new String[] {"forward", "reverse"}) {
              message(" " + mode + " Jacobian on SX. Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
              org.casadi.Function f = new org.casadi.Function("f", List.of(sxinputs.get(inputshape + " " + inputtype)),
                  List.of(sxoutputs.get(outputshape + " " + outputtype)), adOpts(mode));
              org.casadi.Function Jf = jacobian_old(f, 0, 0);
              DM J_in = new DM(f.sparsity_in(0), new DM(n));
              DM Jout = (DM) ((List<?>) Jf.invoke(J_in)).get(0);
              DM J = jacobians(inputtype, outputtype, n);
              checkarray(Jout, J, "Jacobian\n Mode: " + mode + "\n Input: " + inputshape + " " + inputtype + "\n Output: " + outputshape + " " + outputtype);
            }
  }

  @Test
  public void test_jacobianSX() {
    double[] n = {1.2, 2.3, 7, 4.6};
    for (String inputshape : SHAPES)
      for (String outputshape : SHAPES)
        for (String inputtype : TYPES)
          for (String outputtype : TYPES) {
            message("jacobian on SX (SCT). Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
            SX in = sxinputs.get(inputshape + " " + inputtype);
            org.casadi.Function Jf = new org.casadi.Function("Jf", List.of(in),
                List.of(ca.jacobian(new SX(sxoutputs.get(outputshape + " " + outputtype)), new SX(in))));
            DM J_in = new DM(Jf.sparsity_in(0), new DM(n));
            DM J_out = (DM) Jf.invoke(J_in);
            checkarray(J_out, jacobians(inputtype, outputtype, n), "jacobian");
          }
  }

  @Test
  public void test_jacsparsity() {
    double[] n = {1.2, 2.3, 7, 4.6};
    for (String inputshape : SHAPES)
      for (String outputshape : SHAPES)
        for (String inputtype : TYPES)
          for (String outputtype : TYPES) {
            message("jacsparsity on SX. Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
            org.casadi.Function f = new org.casadi.Function("f", List.of(sxinputs.get(inputshape + " " + inputtype)),
                List.of(sxoutputs.get(outputshape + " " + outputtype)));
            DM J = jacobians(inputtype, outputtype, n);
            checkarray(DM.ones(f.jac_sparsity(0, 0)), nonzeroIndicator(J), "jacsparsity");
          }
  }

  @Test
  public void test_JacobianMX() {
    double[] n = {1.2, 2.3, 7, 4.6};
    for (String inputshape : SHAPES)
      for (String outputshape : SHAPES)
        for (String inputtype : TYPES)
          for (String outputtype : TYPES)
            for (String mode : new String[] {"forward", "reverse"}) {
              message("adj AD on MX. Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
              MX in = mxinputs.get(inputshape + " " + inputtype);
              org.casadi.Function f = new org.casadi.Function("f", List.of(in),
                  List.of(mxoutputs.get(outputshape + " " + outputtype).apply(in)), adOpts(mode));
              org.casadi.Function Jf = jacobian_old(f, 0, 0);
              DM J_in = new DM(f.sparsity_in(0), new DM(n));
              DM J_out = (DM) ((List<?>) Jf.invoke(J_in)).get(0);
              DM J = jacobians(inputtype, outputtype, n);
              checkarray(J_out, J, "Jacobian\n Mode: " + mode + "\n Input: " + inputshape + " " + inputtype + "\n Output: " + outputshape + " " + outputtype);
            }
  }

  @Test
  public void test_jacsparsityMX() {
    double[] n = {1.2, 2.3, 7, 4.6};
    for (String inputshape : SHAPES)
      for (String outputshape : SHAPES)
        for (String inputtype : TYPES)
          for (String outputtype : TYPES)
            for (String mode : new String[] {"forward", "reverse"}) {
              message(" " + mode + " jacobian on MX (SCT). Input " + inputtype + " " + inputshape + ", Output " + outputtype + " " + outputshape);
              MX in = mxinputs.get(inputshape + " " + inputtype);
              org.casadi.Function f = new org.casadi.Function("f", List.of(in),
                  List.of(mxoutputs.get(outputshape + " " + outputtype).apply(in)), adOpts(mode));
              org.casadi.Function Jf = jacobian_old(f, 0, 0);
              DM J_in = new DM(f.sparsity_in(0), new DM(n));
              DM J_out = (DM) ((List<?>) Jf.invoke(J_in)).get(0);
              DM J = jacobians(inputtype, outputtype, n);
              checkarray(J_out, J, "jacobian");
              checkarray(DM.ones(f.jac_sparsity(0, 0)), nonzeroIndicator(J), "jacsparsity");
            }
  }

  @Test
  public void test_hessian() {
    message("Jacobian chaining");
    SX x = SX.sym("x"), y = SX.sym("y"), z = SX.sym("z");
    double[] n = {1.2, 2.3, 7};
    org.casadi.Function f = new org.casadi.Function("f", List.of(ca.vertcat(List.of(x, y, z))),
        List.of(ca.vertcat(List.of(x.plus(ca.times(2.0, y.power(3.0))).plus(ca.times(3.0, z.power(4.0)))))));
    org.casadi.Function J = jacobian_old(f, 0, 0);
    MX m = MX.sym("m", 3, 1);
    MX JT = (MX) ((List<?>) J.invoke(m)).get(0);
    org.casadi.Function JTf = new org.casadi.Function("JT", List.of(m), List.of(JT.T()));
    JTf.invoke(n);
    org.casadi.Function H = jacobian_old(JTf, 0, 0);
    H.invoke(n);
  }

  @Test
  public void test_bugshape() {
    message("shape bug");
    SX x = SX.sym("x"), y = SX.sym("y");
    SX inp = new SX(5, 1);
    inp.set(0, 0, x);
    inp.set(3, 0, y);
    org.casadi.Function f = new org.casadi.Function("f", List.of(inp), List.of(ca.vertcat(List.of(x.plus(y), x, y))));
    org.casadi.Function J = jacobian_old(f, 0, 0);
    J.invoke(new DM(f.sparsity_in(0), dm(2, 7)));
    assertEquals(3, f.size1_out(0), "Jacobian shape bug");
    assertEquals(1, f.size2_out(0), "Jacobian shape bug");
  }

  @Test
  public void test_bugglibc() {
    message("Code that used to throw a glibc error");
    SX x = SX.sym("x"), y = SX.sym("y");
    SX inp = new SX(5, 1);
    inp.set(0, 0, x);
    inp.set(3, 0, y);
    org.casadi.Function f = new org.casadi.Function("f", List.of(inp), List.of(ca.vertcat(List.of(x.plus(y), x, y))));
    org.casadi.Function J = jacobian_old(f, 0, 0);
    DM J_in = new DM(f.sparsity_in(0), dm(2, 7));
    ((List<?>) J.invoke(J_in)).get(0);
    f = new org.casadi.Function("f", List.of(inp), List.of(ca.vertcat(List.of(x.plus(y), x, y))));
    jacobian_old(f, 0, 0);
  }

  private static MX mx(Object o) {
    if (o instanceof MX) return (MX) o;
    if (o instanceof DM) return new MX((DM) o);
    return new MX(((Number) o).doubleValue());
  }

  private static MX bc(Object[]... rows) {
    List<List<MX>> r = new ArrayList<>();
    for (Object[] row : rows) {
      List<MX> l = new ArrayList<>();
      for (Object o : row) l.add(mx(o));
      r.add(l);
    }
    return ca.blockcat(r);
  }

  private static DM idx(double... i) {
    return new DM(i);
  }

  private static MX sq(MX a) {
    return a.power(2.0);
  }

  private static final class Case {
    final List<DM> values;
    final MX out, jac;
    final boolean withSx;
    final String std;

    Case(List<DM> values, Object out, Object jac, boolean withSx, String std) {
      this.values = values;
      this.out = mx(out);
      this.jac = mx(jac);
      this.withSx = withSx;
      this.std = std;
    }
  }

  private static List<DM> vecl(List<List<DM>> l) {
    List<DM> r = new ArrayList<>();
    for (List<DM> i : l) r.addAll(i);
    return r;
  }

  private static <T> List<T> concat(List<? extends T> a, List<? extends T> b) {
    List<T> r = new ArrayList<>(a);
    r.addAll(b);
    return r;
  }

  private static DM randDMsp(Random rnd, Sparsity sp) {
    double[] v = new double[(int) sp.nnz()];
    for (int k = 0; k < v.length; k++) v[k] = rnd.nextDouble();
    return new DM(sp, new DM(v));
  }

  private Map<String, List<Object[]>> cgBatch;

  private void check_codegen(org.casadi.Function F, List<DM> inputs, String std) {
    if (!RUN_SLOW) return;
    cgBatch.computeIfAbsent(std, k -> new ArrayList<>()).add(new Object[] {F, inputs});
  }

  private static int cgCounter = 0;

  private void flushCodegen() throws Exception {
    Map<String, List<Object[]>> batch = cgBatch;
    cgBatch = new LinkedHashMap<>();
    for (Map.Entry<String, List<Object[]>> e : batch.entrySet()) {
      String name = "codegen_ad_" + ProcessHandle.current().pid() + "_" + (cgCounter++);
      File dir = Files.createTempDirectory("casadi_ad").toFile();
      CodeGenerator cg = new CodeGenerator(name);
      List<org.casadi.Function> wrapped = new ArrayList<>();
      for (int k = 0; k < e.getValue().size(); k++) {
        org.casadi.Function f = (org.casadi.Function) e.getValue().get(k)[0];
        org.casadi.Function fw = f.wrap(f.name() + "_b" + k);
        wrapped.add(fw);
        cg.add(fw);
      }
      cg.generate(dir.getPath() + "/");
      String inc = GlobalOptions.getCasadiIncludePath(), lib = GlobalOptions.getCasadiPath();
      String so = dir.getPath() + "/" + name + ".so";
      Process p = new ProcessBuilder("gcc", "-pedantic", "-std=" + e.getKey(), "-fPIC", "-shared", "-Wall", "-Werror", "-Wextra",
          "-I" + inc, "-Wno-dangling-pointer", "-Wno-unknown-warning-option", "-Wno-unknown-pragmas", "-Wno-long-long",
          "-Wno-unused-parameter", "-Wno-misleading-indentation", "-O3", dir.getPath() + "/" + name + ".c", "-o", so,
          "-L" + lib, "-Wl,-rpath," + lib).inheritIO().start();
      assertEquals(0, p.waitFor(), "gcc failed");
      for (int k = 0; k < wrapped.size(); k++) {
        @SuppressWarnings("unchecked")
        List<DM> in = (List<DM>) e.getValue().get(k)[1];
        org.casadi.Function f = wrapped.get(k), f2 = ca.external(f.name(), so);
        List<DM> r = f.call(in), r2 = f2.call(in);
        for (int i = 0; i < r.size(); i++) checkarray(r.get(i), r2.get(i), CODEGEN_TOL, "codegen " + f.name());
      }
      for (File g : dir.listFiles()) g.delete();
      dir.delete();
    }
  }

  private void check_serialize(org.casadi.Function F, List<DM> inputs) {
    String s = F.serialize(Map.of("debug", true));
    org.casadi.Function F2 = org.casadi.Function.deserialize(s);
    List<DM> a = F.call(inputs), b = F2.call(inputs);
    for (int i = 0; i < F.n_out(); i++) checkarray(a.get(i), b.get(i), 0.5e-16, "serialize");
  }

  private void check_eval_mx(List<MX> l) {
    check_eval_mx(ca.veccat(l));
  }

  private void check_eval_mx(MX mx) {
    List<MX> args = ca.symvar(mx);
    org.casadi.Function f = new org.casadi.Function("f", args, List.of(mx));
    List<MX> new_args = new ArrayList<>();
    for (MX e : args) new_args.add(MX.sym(e.name(), e.sparsity()));
    Object mx_eval = ((List<?>) f.call(new_args, true, false)).get(0);
    if (mx_eval instanceof DM) {
      checkarray(ca.evalf(mx), (DM) mx_eval);
      return;
    }
    MX ev = ca.substitute(List.of((MX) mx_eval), new_args, args).get(0);
    check_identical_mx(ev, mx);
  }

  private static void check_identical_mx(MX a, MX b) {
    org.casadi.Function fa = new org.casadi.Function("f", ca.symvar(a), List.of(a));
    org.casadi.Function fb = new org.casadi.Function("f", ca.symvar(b), List.of(b));
    assertTrue(fa.serialize().equals(fb.serialize()), "identical mx");
  }

  private static Sparsity remove_first(Sparsity x) {
    DM ret = new DM(x);
    if (ret.numel() > 0) ret.set(0, 0, new DM(1, 1));
    return ret.sparsity();
  }

  private static Sparsity remove_last(Sparsity x) {
    DM ret = new DM(x);
    if (ret.nnz() > 0) {
      long[] r = ret.sparsity().row(), c = ret.sparsity().get_col();
      ret.set(r[r.length - 1], c[c.length - 1], new DM(1, 1));
      return ret.sparsity();
    }
    return x;
  }

  @SuppressWarnings("unchecked")
  @Test
  public void test_MX() throws Exception {
    MX x = MX.sym("x", 2), y = MX.sym("y", 2, 2);
    org.casadi.Function f1 = new org.casadi.Function("f1", List.of(x, y), List.of(x.plus(y.get(0, 0)), ca.mtimes(y, x)));
    org.casadi.Function f1_noninline = new org.casadi.Function("f1", List.of(x, y), List.of(x.plus(y.get(0, 0)), ca.mtimes(y, x)),
        Map.of("never_inline", true));
    org.casadi.Function f2 = new org.casadi.Function("f2", List.of(x, y), List.of(ca.mtimes(MX.zeros(0, 2), x)));
    org.casadi.Function f3 = new org.casadi.Function("f3", List.of(x, y), List.of(MX.zeros(0, 0), ca.mtimes(y, x)));
    org.casadi.Function f4 = new org.casadi.Function("f4", List.of(x, y), List.of(MX.zeros(0, 2), ca.mtimes(y, x)));
    int ndir = 2;
    List<MX> in1 = List.of(x, y);
    List<DM> v1 = List.of(dm(1.1, 1.3), dm(new double[][] {{0.7, 1.5}, {2.1, 0.9}}));
    Slice all = new Slice();
    MX w = x.get(all);
    w.set(1, w.get(1).times(2.0));
    MX w2 = x.get(all);
    w2.set(1, w2.get(1).times(x.get(0)));
    MX ww = x.get(all);
    ww.set(idx(0, 1), ww.get(idx(0, 1)).times(x));
    MX wwf = x.get(all);
    wwf.set(idx(1, 0), wwf.get(idx(1, 0)).times(x));
    MX wwr = x.get(all);
    wwr.set(idx(0, 0, 1, 1), wwr.get(idx(0, 0, 1, 1)).times(2.0));
    MX yy = y.get(all, all);
    yy.set(all, 0, x);
    MX yy2 = y.get(all, all);
    yy2.set(all, 0, sq(x));
    MX yyy = y.get(all, all);
    yyy.set(idx(1, 0), 0, x);
    MX yyy2 = y.get(all, all);
    yyy2.set(idx(1, 0), 0, sq(x));
    List<UnaryOperator<Sparsity>> spmods = List.of(s -> s, ADTest::remove_first, ADTest::remove_last);
    MX y_nz = ca.sin(y);
    y_nz.set_nz(ca.vertcat(List.of(x.get(1), x.get(0))), false, idx(2, 1));
    MX y_nzx = ca.sin(x);
    y_nzx.set_nz(sq(x.get(0)), false, 1);
    MX y_nzs = ca.sin(y);
    y_nzs.set_nz(ca.vertcat(List.of(x.get(1), x.get(0))), false, y.get(all, 0));
    MX y_nzxs = ca.sin(x);
    y_nzxs.set_nz(sq(x.get(0)), false, y.get(1, 0));
    MX xx = ca.horzcat(List.of(ca.sin(x), ca.cos(x)));
    MX Z = new MX(1, 1);
    DM eye2 = DM.eye(2);
    List<DM> v2 = List.of(v1.get(0), dm(new double[][] {{1, 1.5}, {0, 0.9}}));
    Sparsity t1 = Sparsity.triplet(2, 1, new long[] {1}, new long[] {0});
    Sparsity t2 = Sparsity.triplet(2, 2, new long[] {1, 1}, new long[] {0, 1});
    Sparsity t3 = Sparsity.triplet(2, 2, new long[] {0, 1, 1}, new long[] {0, 0, 1});
    Object[][] blk16 = new Object[16][];
    blk16[0] = new Object[] {Z, x.get(1).times(2.0)};
    for (int i = 1; i < 15; i++) blk16[i] = new Object[] {Z, Z};
    blk16[15] = new Object[] {x.get(0).times(2.0), Z};
    List<Case> cases = List.of(
        new Case(v1, ca.sparsity_cast(sq(x), ca.sparsify(dm(new double[][] {{0, 1}, {1, 0}})).sparsity()),
            ca.sparsity_cast(x, ca.sparsify(dm(new double[][] {{0, 0}, {1, 0}, {0, 1}, {0, 0}})).sparsity()).times(2.0), true, "c89"),
        new Case(v1, f1_noninline.call(List.of(sq(x), y)).get(1), y.times(2.0).times(ca.vertcat(List.of(x.T(), x.T()))), true, "c89"),
        new Case(v2, xx.get(false, y.get(all, 0), all), bc(new Object[] {0, ca.cos(x.get(1))}, new Object[] {ca.cos(x.get(0)), 0},
            new Object[] {0, ca.sin(x.get(1)).unaryMinus()}, new Object[] {ca.sin(x.get(0)).unaryMinus(), 0}), false, "c99"),
        new Case(v2, xx.get(false, all, y.get(all, 0)), bc(new Object[] {ca.sin(x.get(0)).unaryMinus(), 0},
            new Object[] {0, ca.sin(x.get(1)).unaryMinus()}, new Object[] {ca.cos(x.get(0)), 0}, new Object[] {0, ca.cos(x.get(1))}), false, "c99"),
        new Case(v2, xx.get(false, y.get(all, 0), y.get(all, 0)), bc(new Object[] {0, ca.sin(x.get(1)).unaryMinus()},
            new Object[] {ca.sin(x.get(0)).unaryMinus(), 0}, new Object[] {0, ca.cos(x.get(1))}, new Object[] {ca.cos(x.get(0)), 0}), false, "c99"),
        new Case(v1, y_nz, bc(new Object[] {Z, Z}, new Object[] {1, Z}, new Object[] {Z, 1}, new Object[] {Z, Z}), true, "c89"),
        new Case(v1, y_nzx, bc(new Object[] {ca.cos(x.get(0)), Z}, new Object[] {x.get(0).times(2.0), Z}), true, "c89"),
        new Case(List.of(dm(2, 1), v1.get(1)), y.get_nz(false, x), new DM(2, 2), false, "c99"),
        new Case(v2, ca.sin(x).get_nz(false, y.get(all, 0)), bc(new Object[] {0, ca.cos(x.get(1))}, new Object[] {ca.cos(x.get(0)), 0}), false, "c99"),
        new Case(List.of(v1.get(0), dm(new double[][] {{2, 1.5}, {1, 0.9}})), y_nzs,
            bc(new Object[] {0, 0}, new Object[] {1, 0}, new Object[] {0, 1}, new Object[] {0, 0}), false, "c99"),
        new Case(List.of(v1.get(0), dm(new double[][] {{0, 1.5}, {1, 0.9}})), y_nzxs,
            bc(new Object[] {ca.cos(x.get(0)), 0}, new Object[] {x.get(0).times(2.0), 0}), false, "c99"),
        new Case(List.of(v1.get(0), dm(new double[][] {{1, 1.5}, {1, 0.9}})), x.get_nz(false, y.get(all, 0)),
            bc(new Object[] {0, 1}, new Object[] {0, 1}), false, "c99"),
        new Case(v1, x, eye2, true, "c89"),
        new Case(v1, x.T(), eye2, true, "c89"),
        new Case(v1, sq(x), ca.diag(x).times(2.0), true, "c89"),
        new Case(v1, sq(x).attachAssert(1.0), ca.diag(x).times(2.0), true, "c89"),
        new Case(v1, sq(x).T(), ca.diag(x).times(2.0), true, "c89"),
        new Case(v1, ca.reshape(x, 1, 2), eye2, true, "c89"),
        new Case(v1, ca.reshape(sq(x), 1, 2), ca.diag(x).times(2.0), true, "c89"),
        new Case(v1, x.plus(y.get_nz(false, 0)), eye2, true, "c89"),
        new Case(v1, x.plus(y.get(0, 0)), eye2, true, "c89"),
        new Case(v1, x.plus(x), eye2.times(2.0), true, "c89"),
        new Case(v1, sq(x).plus(x), ca.diag(x).times(2.0).plus(eye2), true, "c89"),
        new Case(v1, x.times(x), ca.diag(x).times(2.0), true, "c89"),
        new Case(v1, x.times(y.get_nz(false, 0)), ca.times(eye2, y.get_nz(false, 0)), true, "c89"),
        new Case(v1, x.times(y.get(0, 0)), ca.times(eye2, y.get(0, 0)), true, "c89"),
        new Case(v1, x.get(0), eye2.get(0, all), true, "c89"),
        new Case(v1, sq(x).get(0), ca.horzcat(List.of(x.get(0).times(2.0), Z)), true, "c89"),
        new Case(v1, x.get(0).plus(x.get(1)), DM.ones(1, 2), true, "c89"),
        new Case(v1, ca.sin(ca.repmat(sq(x), 1, 3)), ca.repmat(ca.cos(ca.diag(sq(x))).times(2.0).times(ca.diag(x)), 3, 1), true, "c89"),
        new Case(v1, ca.sin(ca.repsum(sq(x).T(), 1, 2)), ca.cos(sq(x.get(0)).plus(sq(x.get(1)))).times(2.0).times(x.T()), true, "c89"),
        new Case(v1, ca.vertcat(List.of(x.get(1), x.get(0))), ca.sparsify(dm(new double[][] {{0, 1}, {1, 0}})), true, "c89"),
        new Case(v1, ca.vertsplit(x, new long[] {0, 1, 2}).get(1), ca.sparsify(dm(new double[][] {{0, 1}})), true, "c89"),
        new Case(v1, ca.vertcat(List.of(sq(x.get(1)), sq(x.get(0)))),
            bc(new Object[] {Z, x.get(1).times(2.0)}, new Object[] {x.get(0).times(2.0), Z}), true, "c89"),
        new Case(v1, ca.vertsplit(sq(x), new long[] {0, 1, 2}).get(1), bc(new Object[] {Z, x.get(1).times(2.0)}), true, "c89"),
        new Case(v1, ca.vertsplit(sq(x), new long[] {0, 1, 2}).get(1).power(3.0), bc(new Object[] {Z, x.get(1).power(5.0).times(6.0)}), true, "c89"),
        new Case(v1, ca.horzcat(List.of(x.get(1), x.get(0))).T(), ca.sparsify(dm(new double[][] {{0, 1}, {1, 0}})), true, "c89"),
        new Case(v1, ca.horzcat(List.of(sq(x.get(1)), sq(x.get(0)))).T(),
            bc(new Object[] {Z, x.get(1).times(2.0)}, new Object[] {x.get(0).times(2.0), Z}), true, "c89"),
        new Case(v1, ca.diagcat(List.of(sq(x.get(1)), y, sq(x.get(0)))), bc(blk16), true, "c89"),
        new Case(v1, ca.horzcat(List.of(sq(x.get(1)), sq(x.get(0)))).T(),
            bc(new Object[] {Z, x.get(1).times(2.0)}, new Object[] {x.get(0).times(2.0), Z}), true, "c89"),
        new Case(v1, x.get(idx(0, 1)), ca.sparsify(dm(new double[][] {{1, 0}, {0, 1}})), true, "c89"),
        new Case(v1, sq(x).get(idx(0, 1)), ca.diag(x).times(2.0), true, "c89"),
        new Case(v1, x.get(idx(0, 0, 1, 1)), ca.sparsify(dm(new double[][] {{1, 0}, {1, 0}, {0, 1}, {0, 1}})), true, "c89"),
        new Case(v1, sq(x).get(idx(0, 0, 1, 1)), bc(new Object[] {x.get(0).times(2.0), Z}, new Object[] {x.get(0).times(2.0), Z},
            new Object[] {Z, x.get(1).times(2.0)}, new Object[] {Z, x.get(1).times(2.0)}), true, "c89"),
        new Case(v1, wwr, ca.sparsify(dm(new double[][] {{2, 0}, {0, 2}})), true, "c89"),
        new Case(v1, x.get(idx(1, 0)), ca.sparsify(dm(new double[][] {{0, 1}, {1, 0}})), true, "c89"),
        new Case(v1, x.get(idx(1, 0), 0), ca.sparsify(dm(new double[][] {{0, 1}, {1, 0}})), true, "c89"),
        new Case(v1, w, ca.sparsify(dm(new double[][] {{1, 0}, {0, 2}})), true, "c89"),
        new Case(v1, w2, bc(new Object[] {1, Z}, new Object[] {x.get(1), x.get(0)}), true, "c89"),
        new Case(v1, ww, ca.diag(x).times(2.0), true, "c89"),
        new Case(v1, wwf, ca.vertcat(List.of(x.get(idx(1, 0)).T(), x.get(idx(1, 0)).T())), true, "c89"),
        new Case(v1, yy.get(all, 0), eye2, true, "c89"),
        new Case(v1, yy2.get(all, 0), ca.diag(x).times(2.0), true, "c89"),
        new Case(v1, yyy.get(all, 0), ca.sparsify(dm(new double[][] {{0, 1}, {1, 0}})), true, "c89"),
        new Case(v1, ca.mtimes(y, x), y, true, "c89"),
        new Case(v1, ca.mtimes(x.T(), y.T()), y, true, "c89"),
        new Case(v1, ca.mac(y, x, DM.zeros(t1)), y.get(t2), true, "c89"),
        new Case(v1, ca.mac(x.T(), y.T(), DM.zeros(t1.T())), y.get(t2), true, "c89"),
        new Case(v1, ca.mtimes(y.get(t3), x), y.get(t3), true, "c89"),
        new Case(v1, ca.mtimes(x.T(), y.get(t3).T()), y.get(t3), true, "c89"),
        new Case(v1, ca.mtimes(y, sq(x)), y.times(2.0).times(ca.vertcat(List.of(x.T(), x.T()))), true, "c89"),
        new Case(v1, ca.sin(x), ca.diag(ca.cos(x)), true, "c89"),
        new Case(v1, ca.sin(sq(x)), ca.diag(ca.cos(sq(x)).times(2.0).times(x)), true, "c89"),
        new Case(v1, x.times(y.get(all, 0)), ca.diag(y.get(all, 0)), true, "c89"),
        new Case(v1, x.times(y.get_nz(false, idx(0, 1))), ca.diag(y.get_nz(false, idx(0, 1))), true, "c89"),
        new Case(v1, x.times(y.get_nz(false, idx(1, 0))), ca.diag(y.get_nz(false, idx(1, 0))), true, "c89"),
        new Case(v1, x.times(y.get(idx(0, 1), 0)), ca.diag(y.get(idx(0, 1), 0)), true, "c89"),
        new Case(v1, x.times(y.get(idx(1, 0), 0)), ca.diag(y.get(idx(1, 0), 0)), true, "c89"),
        new Case(v1, ca.dot(x, x), x.times(2.0).T(), true, "c89"),
        new Case(v1, ca.dot(sq(x), x), sq(x).times(3.0).T(), true, "c89"),
        new Case(v1, ca.det(ca.horzcat(List.of(x, new MX(dm(1, 2))))), dm(new double[][] {{2, -1}}), true, "c89"),
        new Case(v1, ca.det(ca.horzcat(List.of(ca.vertcat(List.of(x.get(0), new MX(1.0))), ca.vertcat(List.of(new MX(1.0), x.get(1)))))),
            ca.horzcat(List.of(x.get(1), x.get(0))), true, "c89"),
        new Case(v1, f1.call(in1).get(1), y, true, "c89"),
        new Case(v1, f1.call(List.of(sq(x), y)).get(1), y.times(2.0).times(ca.vertcat(List.of(x.T(), x.T()))), true, "c89"),
        new Case(v1, f2.call(in1).get(0), DM.zeros(0, 2), true, "c89"),
        new Case(v1, f2.invoke(sq(x), y), DM.zeros(0, 2), true, "c89"),
        new Case(v1, f3.call(in1).get(0), DM.zeros(0, 2), true, "c89"),
        new Case(v1, f3.call(List.of(sq(x), y)).get(0), DM.zeros(0, 2), true, "c89"),
        new Case(v1, f4.call(in1).get(0), DM.zeros(0, 2), true, "c89"),
        new Case(v1, f4.call(List.of(sq(x), y)).get(0), DM.zeros(0, 2), true, "c89"),
        new Case(v1, ca.vertcat(List.of(x, new MX(new DM(0, 1)))), eye2, true, "c89"),
        new Case(v1, ca.project(sq(x), ca.sparsify(dm(0, 1)).sparsity()), bc(new Object[] {Z, Z}, new Object[] {Z, x.get(1).times(2.0)}), true, "c89"),
        new Case(v1, ca.dot(x, y.get(all, 0)), y.get(all, 0).T(), true, "c89"),
        new Case(v1, x.get_nz(false, dm(new double[][] {{1, 0}})).T().times(y.get_nz(false, dm(new double[][] {{0, 2}}))),
            bc(new Object[] {Z, y.get_nz(false, 0)}, new Object[] {y.get_nz(false, 2), Z}), true, "c89"),
        new Case(v1, x.get_nz(false, ca.diag(dm(1, 0))).times(y.get_nz(false, ca.diag(dm(0, 2)))),
            bc(new Object[] {Z, y.get_nz(false, 0)}, new Object[] {Z, Z}, new Object[] {Z, Z}, new Object[] {y.get_nz(false, 2), Z}), true, "c89"),
        new Case(v1, ca.kron(x, y), ca.horzcat(List.of(ca.vec(ca.vertcat(List.of(y, new MX(2, 2)))), ca.vec(ca.vertcat(List.of(new MX(2, 2), y))))), true, "c89"),
        new Case(v1, ca.kron(y, x), ca.vertcat(List.of(ca.times(y.get(0, 0), eye2), ca.times(y.get(1, 0), eye2),
            ca.times(y.get(0, 1), eye2), ca.times(y.get(1, 1), eye2))), true, "c89"),
        new Case(v1, ca.kron_contract(ca.kron(x, y), y, true), ca.times(ca.sumsqr(y), eye2), true, "c89"),
        new Case(v1, ca.kron_contract(ca.kron(x, y), x, false), ca.mtimes(ca.vec(y).times(2.0), x.T()), true, "c89"));

    cgBatch = new LinkedHashMap<>();
    for (Case cs : cases) {
      List<DM> values = cs.values;
      MX jac = cs.jac;
      boolean with_sx = cs.withSx;
      String std = cs.std;
      org.casadi.Function fun0 = new org.casadi.Function("fun", in1, List.of(cs.out));
      List<MX> shifted = new ArrayList<>();
      for (MX e : in1) shifted.add(e.plus(1e-300));
      MX out_eval_mx = fun0.call(shifted, true, false).get(0);
      for (MX out : List.of(cs.out, out_eval_mx)) {
        check_eval_mx(out);
        org.casadi.Function fun = new org.casadi.Function("fun", in1, List.of(out, jac));
        List<org.casadi.Function> fun_ad = new ArrayList<>();
        for (int wt : new int[] {0, 1})
          fun_ad.add(new org.casadi.Function("fun", in1, List.of(out, jac), Map.of("ad_weight", wt, "ad_weight_sp", wt)));
        List<DM> fun_out = fun.call(values);
        List<org.casadi.Function> funsx_ad = null;
        if (with_sx) {
          org.casadi.Function funsx = fun.expand("expand_fun");
          funsx_ad = new ArrayList<>();
          for (org.casadi.Function f : fun_ad) funsx_ad.add(f.expand("expand_" + f.name()));
          List<DM> funsx_out = funsx.call(values);
          checkarray(fun_out.get(0), funsx_out.get(0));
          checkarray(fun_out.get(1), funsx_out.get(1));
        }
        check_codegen(fun, values, std);
        check_serialize(fun, values);
        DM J_ = fun_out.get(1);
        Map<String, List<List<DM>>> storage2 = new LinkedHashMap<>(), storage = new LinkedHashMap<>();
        org.casadi.Function vf_mx = null;
        List<org.casadi.Function> fs = new ArrayList<>(List.of(fun));
        if (with_sx) fs.add(fun.expand("expand_" + fun.name()));
        for (org.casadi.Function f : fs) {
          f.forward(ndir);
          f.reverse(ndir);
          for (boolean symMX : new boolean[] {true, false}) {
            if (f.is_a("MXFunction") && !symMX) continue;
            if (f.is_a("SXFunction") && symMX) continue;
            BiFunction<String, Sparsity, GenericMatrix> sym = symMX ? MX::sym : SX::sym;
            org.casadi.Function vf = null;
            String storagekey = null;
            for (int a = 0; a < 3; a++)
              for (int b = 0; b < 3; b++) {
                UnaryOperator<Sparsity> spmod = spmods.get(a), spmod2 = spmods.get(b);
                storagekey = a + "," + b;
                List<List<GenericMatrix>> fseeds = new ArrayList<>(), aseeds = new ArrayList<>();
                for (int d = 0; d < ndir; d++) {
                  List<GenericMatrix> l = new ArrayList<>();
                  for (long i = 0; i < f.n_in(); i++) l.add(sym.apply("f", spmod.apply(f.sparsity_in(i))));
                  fseeds.add(l);
                  l = new ArrayList<>();
                  for (long i = 0; i < f.n_out(); i++) l.add(sym.apply("a", spmod2.apply(f.sparsity_out(i))));
                  aseeds.add(l);
                }
                List<GenericMatrix> inputss = new ArrayList<>();
                for (long i = 0; i < f.n_in(); i++) inputss.add(sym.apply("i", f.sparsity_in(i)));
                boolean notSX = !f.is_a("SXFunction");
                List<GenericMatrix> res = f.call(inputss, notSX);
                List<List<GenericMatrix>> fwdsens = ca.forward(res, inputss, fseeds, Map.of("always_inline", notSX));
                List<List<GenericMatrix>> adjsens = ca.reverse(res, inputss, aseeds, Map.of("always_inline", notSX));
                List<DM> fseed = new ArrayList<>(), aseed = new ArrayList<>();
                for (int d = 0; d < ndir; d++) {
                  fseed.add(randDMsp(random, spmod.apply(f.sparsity_in(0))));
                  aseed.add(randDMsp(random, spmod2.apply(f.sparsity_out(0))));
                }
                List<GenericMatrix> vin = new ArrayList<>(inputss), vout = new ArrayList<>(res);
                for (int d = 0; d < ndir; d++) {
                  vin.addAll(fseeds.get(d));
                  vin.addAll(aseeds.get(d));
                  vout.addAll(fwdsens.get(d));
                  vout.addAll(adjsens.get(d));
                }
                vf = new org.casadi.Function("vf", vin, vout);
                List<DM> vf_in = new ArrayList<>(values);
                for (int d = 0; d < ndir; d++) {
                  vf_in.add(fseed.get(d));
                  for (int i = 0; i < values.size() - 1; i++) vf_in.add(new DM(0.0));
                  vf_in.add(aseed.get(d));
                  vf_in.add(new DM(0.0));
                }
                List<DM> vf_out = vf.call(vf_in);
                if (symMX) {
                  List<MX> ef = new ArrayList<>(), ea = new ArrayList<>();
                  for (List<GenericMatrix> e : fwdsens) ef.add(ca.veccat((List<MX>) (List<?>) e));
                  for (List<GenericMatrix> e : adjsens) ea.add(ca.veccat((List<MX>) (List<?>) e));
                  check_eval_mx(ef);
                  check_eval_mx(ea);
                }
                check_codegen(vf, vf_in, std);
                check_serialize(vf, vf_in);
                int offset = res.size();
                for (int d = 0; d < ndir; d++) {
                  DM sens = vf_out.get(offset);
                  offset += inputss.size();
                  checkarray(ca.vec(sens), ca.mtimes(J_, ca.vec(fseed.get(d))), "eval Fwd " + d + " " + f.class_name() + " " + symMX);
                  sens = vf_out.get(offset);
                  offset += inputss.size();
                  checkarray(ca.vec(sens), ca.mtimes(J_.T(), ca.vec(aseed.get(d))), "eval Adj " + d + " " + vf_out);
                }
                assertEquals(vf.n_out(), offset);
                Random rnd = new Random(1);
                vf_in = new ArrayList<>();
                for (long i = 0; i < vf.n_in(); i++) vf_in.add(randDMsp(rnd, vf.sparsity_in(i)));
                vf_out = vf.call(vf_in);
                check_codegen(vf, vf_in, std);
                check_serialize(vf, vf_in);
                storage.computeIfAbsent(storagekey, k -> new ArrayList<>()).add(vf_out);
                if (symMX) vf_mx = vf;
              }
            for (boolean sym2MX : new boolean[] {true, false}) {
              if (vf.is_a("MXFunction") && !sym2MX) continue;
              if (vf.is_a("MXFunction") && sym2MX) continue;
              BiFunction<String, Sparsity, GenericMatrix> sym2 = sym2MX ? MX::sym : SX::sym;
              for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++) {
                  List<List<GenericMatrix>> fseeds2 = new ArrayList<>(), aseeds2 = new ArrayList<>();
                  for (int d = 0; d < ndir; d++) {
                    List<GenericMatrix> l = new ArrayList<>();
                    for (long i = 0; i < vf.n_in(); i++) l.add(sym2.apply("f", vf_mx.sparsity_in(i)));
                    fseeds2.add(l);
                    l = new ArrayList<>();
                    for (long i = 0; i < vf.n_out(); i++) l.add(sym2.apply("a", vf_mx.sparsity_out(i)));
                    aseeds2.add(l);
                  }
                  List<GenericMatrix> inputss2 = new ArrayList<>();
                  for (long i = 0; i < vf.n_in(); i++) inputss2.add(sym2.apply("i", vf_mx.sparsity_in(i)));
                  boolean notSX = !vf.is_a("SXFunction");
                  List<GenericMatrix> res2 = vf.call(inputss2, notSX);
                  List<List<GenericMatrix>> fwdsens2 = ca.forward(res2, inputss2, fseeds2, Map.of("always_inline", notSX));
                  List<List<GenericMatrix>> adjsens2 = ca.reverse(res2, inputss2, aseeds2, Map.of("always_inline", notSX));
                  if (sym2MX) {
                    List<MX> ef = new ArrayList<>(), ea = new ArrayList<>();
                    for (List<GenericMatrix> e : fwdsens2) ef.add(ca.veccat((List<MX>) (List<?>) e));
                    for (List<GenericMatrix> e : adjsens2) ea.add(ca.veccat((List<MX>) (List<?>) e));
                    check_eval_mx(ef);
                    check_eval_mx(ea);
                  }
                  List<GenericMatrix> vin = new ArrayList<>(inputss2), vout = new ArrayList<>(res2);
                  for (int d = 0; d < ndir; d++) {
                    vin.addAll(fseeds2.get(d));
                    vin.addAll(aseeds2.get(d));
                    vout.addAll(fwdsens2.get(d));
                    vout.addAll(adjsens2.get(d));
                  }
                  org.casadi.Function vf2 = new org.casadi.Function("vf2", vin, vout);
                  Random rnd = new Random(1);
                  List<DM> vf2_in = new ArrayList<>();
                  for (long i = 0; i < vf2.n_in(); i++) vf2_in.add(randDMsp(rnd, vf2.sparsity_in(i)));
                  List<DM> vf2_out = vf2.call(vf2_in);
                  check_codegen(vf2, vf2_in, std);
                  check_serialize(vf2, vf2_in);
                  storage2.computeIfAbsent(storagekey, k -> new ArrayList<>()).add(vf2_out);
                }
            }
          }
        }
        for (Map<String, List<List<DM>>> store : List.of(storage, storage2)) {
          String order = store == storage ? "first-order" : "second-order";
          for (List<List<DM>> st : store.values())
            for (int i = 0; i < st.size() - 1; i++)
              for (int k = 0; k < st.get(0).size(); k++) {
                DM a = st.get(0).get(k), b = st.get(i + 1).get(k);
                if (b.numel() == 0 && ca.sparsify(a).nnz() == 0) continue;
                if (a.numel() == 0 && ca.sparsify(b).nnz() == 0) continue;
                checkarray(ca.sparsify(a), ca.sparsify(b), order + ", output(" + k + ")");
              }
        }
        for (boolean expand : with_sx ? new boolean[] {false, true} : new boolean[] {false})
          for (String mode : new String[] {"forward", "reverse"}) {
            int ind = mode.equals("forward") ? 0 : 1;
            org.casadi.Function f = expand ? funsx_ad.get(ind) : fun_ad.get(ind);
            org.casadi.Function Jf = jacobian_old(f, 0, 0);
            List<DM> Jf_out = Jf.call(values);
            check_codegen(Jf, values, std);
            check_serialize(Jf, values);
            checkarray(Jf_out.get(0), J_);
            checkarray(DM.ones(Jf.sparsity_out(0)), DM.ones(J_.sparsity()), out + mode);
            checkarray(DM.ones(f.jac_sparsity(0, 0)), DM.ones(J_.sparsity()));
          }
        if (out.is_empty()) continue;
        long s_i = out.sparsity().row()[0], s_j = out.sparsity().get_col()[0];
        long s_k = s_i * out.size2() + s_j;
        DM H_ = null;
        for (boolean expand : with_sx ? new boolean[] {false, true} : new boolean[] {false})
          for (String mode : new String[] {"forward", "reverse"}) {
            int wt = mode.equals("forward") ? 0 : 1;
            org.casadi.Function f = new org.casadi.Function("fun", in1, List.of(out.get(s_i, s_j), jac.get(s_k, all).T()),
                Map.of("ad_weight", wt, "ad_weight_sp", wt));
            if (expand) f = f.expand("expand_" + f.name());
            f.call(values);
            org.casadi.Function Hf = hessian_old(f, 0, 0);
            List<DM> Hf_out = Hf.call(values);
            check_codegen(Hf, values, std);
            check_serialize(Hf, values);
            if (H_ == null) H_ = Hf_out.get(0);
            checkarray(Hf_out.get(0), H_, "mode: " + mode);
          }
      }
      flushCodegen();
    }
  }

  @Test
  public void test_repmat() {
    MX X = MX.sym("x", 2, 2);
    for (MX e : List.of(X.T().plus(ca.repmat(MX.zeros(2, 1), 1, 2))))
      for (int weight : new int[] {0, 1}) {
        org.casadi.Function F = new org.casadi.Function("f", List.of(ca.vec(X)), List.of(e), Map.of("ad_weight_sp", weight));
        Sparsity J0 = F.jac_sparsity(0, 0, false);
        J0.spy();
        Sparsity J1 = F.expand().jac_sparsity(0, 0, false);
        J1.spy();
        assertTrue(J0.is_equal(J1));
      }
  }
}
