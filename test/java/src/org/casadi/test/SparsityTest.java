package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GlobalOptions;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.Sparsity;
import org.casadi.ca;
import org.junit.jupiter.api.Test;

/** Port of test/python/sparsity.py */
public class SparsityTest extends CasadiTestCase {

  private static Sparsity fromNz(long n, long m, long[][] nz) {
    Sparsity a = new Sparsity(n, m);
    for (long[] i : nz) a.add_nz(i[0], i[1]);
    return a;
  }

  private static Set<List<Long>> nzSet(long[][] nz) {
    Set<List<Long>> r = new HashSet<>();
    for (long[] i : nz) r.add(List.of(i[0], i[1]));
    return r;
  }

  private static DM idx(long... i) {
    double[] r = new double[i.length];
    for (int k = 0; k < i.length; k++) r[k] = i[k];
    return new DM(r);
  }

  private static long[] longs(Object o) {
    return (long[]) o;
  }

  /** Python helper randDM(n, m, sparsity, symm=True) */
  private DM randDMSymm(long n, long m, double density) {
    DM d = density < 1 ? randDM(n, m, density) : randDM(n, m);
    return d.plus(d.T()).div(2.0);
  }

  @Test
  public void test_union() {
    message("Sparsity union");
    long[][] nza = {{0, 0}, {0, 1}, {2, 0}, {3, 1}};
    long[][] nzb = {{0, 2}, {0, 0}, {2, 2}};
    Sparsity a = fromNz(4, 5, nza), b = fromNz(4, 5, nzb);
    Sparsity c = a.unite(b);
    c = a.plus(b);
    Set<List<Long>> u = nzSet(nza);
    u.addAll(nzSet(nzb));
    assertEquals(u.size(), c.nnz());
    long[] col = c.get_col();
    for (int k = 0; k < c.nnz(); k++) assertTrue(u.contains(List.of(c.row(k), col[k])));
  }

  @Test
  public void test_intersection() {
    message("Sparsity intersection");
    long[][] nza = {{0, 0}, {0, 1}, {2, 0}, {3, 1}, {2, 3}};
    long[][] nzb = {{0, 2}, {0, 0}, {2, 2}, {2, 3}};
    Sparsity a = fromNz(4, 5, nza), b = fromNz(4, 5, nzb);
    Set<List<Long>> sa = nzSet(nza), sb = nzSet(nzb);
    Sparsity c = a.intersect(b);
    for (int k = 0; k < c.nnz(); k++) {
      List<Long> ind = List.of(c.row(k), c.get_col()[k]);
      assertTrue(sa.contains(ind) && sb.contains(ind));
    }
    c = a.times(b);
    Set<List<Long>> i = new HashSet<>(sa);
    i.retainAll(sb);
    assertEquals(i.size(), c.nnz());
    for (int k = 0; k < c.nnz(); k++) {
      List<Long> ind = List.of(c.row(k), c.get_col()[k]);
      assertTrue(sa.contains(ind) && sb.contains(ind));
    }
  }

  @Test
  public void test_get_nzDense() {
    message("get_nzDense");
    Sparsity a = fromNz(4, 5, new long[][] {{0, 0}, {0, 1}, {2, 0}, {3, 1}});
    DM A = DM.ones(a);
    DM Ad = new DM(A.full());
    for (long i : a.find()) assertEquals(1.0, Ad.get_nz(false, i).to_double());
  }

  @Test
  public void test_find_nonzero() {
    DM d = randDM(20, 10, 0.6);
    Sparsity sp = d.sparsity();
    Sparsity sp2 = Sparsity.nonzeros(20, 10, sp.find());
    assertTrue(sp.is_equal(sp2));
  }

  @Test
  public void test_enlarge() {
    message("enlarge");
    message(":dense");
    Sparsity sp = Sparsity.dense(3, 4);
    long[] col = {1, 2, 4};
    long[] row = {0, 3, 4, 6};
    sp.enlarge(7, 8, col, row);
    double[][] z = new double[7][8];
    for (long i : col)
      for (long j : row) z[(int) i][(int) j] = 1;
    checkarray(DM.ones(sp), z);
    message(":sparse");
    sp = new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1}).T();
    DM n = DM.ones(sp);
    z = new double[7][8];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++) z[(int) col[i]][(int) row[j]] = n.get(i, j).to_double();
    sp.enlarge(7, 8, new long[] {1, 2, 4}, new long[] {0, 3, 4, 6});
    checkarray(DM.ones(sp), z);
  }

  private static DM tomatrix(Sparsity s) {
    DM d = DM.ones(s);
    for (long k = 0; k < d.nnz(); k++) d.set_nz(k + 1, false, k);
    return d;
  }

  @Test
  public void test_NZ() {
    message("NZ constructor");
    long[][] nza = {{0, 0}, {0, 1}, {2, 0}, {2, 3}, {2, 4}, {3, 1}};
    Sparsity a = fromNz(4, 5, nza);
    long[] r = new long[nza.length], c = new long[nza.length];
    for (int k = 0; k < nza.length; k++) {
      r[k] = nza[k][0];
      c[k] = nza[k][1];
    }
    Sparsity b = Sparsity.triplet(4, 5, r, c);
    checkarray(tomatrix(a), tomatrix(b), "rowcol");
  }

  /** First test_rowcol of sparsity.py, shadowed there by the later definition */
  @Test
  public void test_rowcol_shadowed() {
    message("rowcol constructor");
    long[] r = {0, 1, 3};
    long[] c = {1, 4};
    Sparsity a = new Sparsity(4, 5);
    for (long i : r)
      for (long j : c) a.add_nz(i, j);
    Sparsity b = Sparsity.rowcol(r, c, 4, 5);
    checkarray(tomatrix(a), tomatrix(b), "rowcol");
  }

  @Test
  public void test_reshape() {
    message("Reshape");
    Sparsity a = fromNz(4, 5, new long[][] {{0, 0}, {0, 1}, {2, 0}, {2, 3}, {2, 4}, {3, 1}});
    double[][] A = tomatrix(a).full();
    DM B = tomatrix(ca.reshape(a, 2, 10));
    double[][] B_ = new double[2][10];
    for (int k = 0; k < 20; k++) B_[k % 2][k / 2] = A[k % 4][k / 4];
    checkarray(B, B_);
  }

  // test_vec: returns immediately in Python ("doesn't make much sense")

  @Test
  public void test_refcount() {
    DM x = new DM(Sparsity.lower(4), 5);
    Sparsity s = ca.mtimes(x, x).sparsity();
    assertEquals(16, s.numel());
  }

  @Test
  public void test_splower() {
    Sparsity sp = new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1});
    System.out.println(new DM(sp, 1));
    System.out.println(new DM(ca.tril(sp), 1));
    System.out.println(java.util.Arrays.toString(sp.get_lower()));
  }

  @Test
  public void test_diag() {
    message("diag");
    Sparsity A = new Sparsity(5, 5);
    A.add_nz(1, 1);
    A.add_nz(2, 4);
    A.add_nz(3, 3);
    List<Object> r = A.get_diag();
    DM B = DM.ones((Sparsity) r.get(0));
    checkarray(B, new double[][] {{0}, {1}, {0}, {1}, {0}});
    assertArrayEquals(new long[] {0, 1}, longs(r.get(1)));

    A = new Sparsity(5, 1);
    A.add_nz(1, 0);
    A.add_nz(2, 0);
    A.add_nz(4, 0);
    double[][] ref = {{0, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 1}};
    r = A.get_diag();
    checkarray(DM.ones((Sparsity) r.get(0)), ref);
    assertArrayEquals(new long[] {0, 1, 2}, longs(r.get(1)));

    A = new Sparsity(1, 5);
    A.add_nz(0, 1);
    A.add_nz(0, 2);
    A.add_nz(0, 4);
    r = A.get_diag();
    checkarray(DM.ones((Sparsity) r.get(0)), ref);
    assertArrayEquals(new long[] {0, 1, 2}, longs(r.get(1)));
  }

  @Test
  public void test_sparsityindex() {
    message("sparsity indexing");
    long[][] nza = {{0, 0}, {0, 1}, {2, 0}, {2, 3}, {3, 3}, {2, 4}, {3, 1}, {4, 1}};
    Sparsity a = fromNz(5, 5, nza);
    MX b = MX.sym("b", a);
    assertThrows(RuntimeException.class, () -> b.get(Sparsity.diag(3)));
    Sparsity d = Sparsity.diag(5);
    MX c = b.get(d);
    assertTrue(c.sparsity().is_equal(d));
    Function f = new Function("f", List.of(b), List.of(c));
    double[] v = new double[nza.length];
    for (int k = 0; k < v.length; k++) v[k] = k + 1;
    DM fin = new DM(b.sparsity(), new DM(v));
    DM f_out = (DM) f.invoke(fin);
    checkarray(new DM(f_out.nonzeros()), dm(1, 0, 0, 7, 0), "sparsity index");
  }

  @Test
  public void test_get_ccs() {
    message("CCS format");
    Sparsity a = fromNz(4, 5, new long[][] {{0, 0}, {0, 1}, {2, 0}, {2, 3}, {3, 3}, {2, 4}, {3, 1}});
    List<long[]> ab1 = a.get_ccs();
    List<long[]> ab2 = a.T().get_crs();
    System.out.println(java.util.Arrays.toString(ab1.get(0)) + " " + java.util.Arrays.toString(ab1.get(1)));
    System.out.println(java.util.Arrays.toString(ab2.get(0)) + " " + java.util.Arrays.toString(ab2.get(1)));
  }

  private List<Long> randomPermutation(long n) {
    List<Long> p = new ArrayList<>();
    for (long i = 0; i < n; i++) p.add(i);
    Collections.shuffle(p, random);
    return p;
  }

  private static List<Long> blocks(DM Ar) {
    Sparsity ST = Ar.sparsity();
    List<Long> blocks = new ArrayList<>();
    long acc = -1, mc = 0;
    long[] row = ST.row(), colind = ST.colind();
    for (int i = 0; i < Ar.size1(); i++) {
      mc = Math.max(row[(int) colind[i + 1] - 1], mc);
      if (mc == i) {
        blocks.add(i - acc);
        acc = i;
      }
    }
    return blocks;
  }

  @Test
  public void test_dm_diagcat_dense() {
    message("Dulmage-Mendelsohn");
    for (int k = 0; k < 20; k++) {
      List<DM> Ai = new ArrayList<>();
      for (int j = 0; j < 10; j++) {
        long d = 1 + random.nextInt(10);
        Ai.add(randDM(d, d));
      }
      DM A = ca.diagcat(Ai);
      DM perm = idx(toLongs(randomPermutation(A.size1())));
      DM AP = A.get(perm, perm);
      List<Object> r = AP.sparsity().btf();
      DM Ar = AP.get(idx(longs(r.get(1))), idx(longs(r.get(2))));
      List<Long> truth = new ArrayList<>();
      for (DM a : Ai) truth.add(a.size1());
      List<Long> tryme = blocks(Ar);
      Collections.sort(truth);
      Collections.sort(tryme);
      assertEquals(truth, tryme);
    }
  }

  @Test
  public void test_scc_diagcat_sparse() {
    message("scc");
    for (int k = 0; k < 20; k++) {
      List<DM> Ai = new ArrayList<>();
      for (int j = 0; j < 10; j++) {
        long d = 1 + random.nextInt(10);
        Ai.add(randDMSymm(d, d, 0.6));
      }
      DM A = ca.diagcat(Ai);
      DM perm = idx(toLongs(randomPermutation(A.size1())));
      DM AP = A.get(perm, perm);
      List<Object> r = AP.sparsity().scc();
      long n = (Long) r.get(0);
      DM p = idx(longs(r.get(1)));
      DM Ar = AP.get(p, p);
      List<Long> tryme = blocks(Ar);
      assertTrue(n >= Ai.size());
      assertTrue(n >= tryme.size());
    }
  }

  private static void checkBtf(DM A, long ret, long[] rowperm, long[] colperm, long[] rowblock, long[] colblock,
      long[] coarse_rowblock, long[] coarse_colblock) {
    List<Object> r = A.sparsity().btf();
    assertEquals(ret, (long) (Long) r.get(0));
    assertArrayEquals(rowperm, longs(r.get(1)));
    assertArrayEquals(colperm, longs(r.get(2)));
    assertArrayEquals(rowblock, longs(r.get(3)));
    assertArrayEquals(colblock, longs(r.get(4)));
    assertArrayEquals(coarse_rowblock, longs(r.get(5)));
    assertArrayEquals(coarse_colblock, longs(r.get(6)));
  }

  @Test
  public void test_dm() {
    DM A = new DM(6, 4);
    A.set(0, 0, 1);
    A.set(1, 2, 1);
    A.set(2, 2, 1);
    A.set(5, 3, 1);
    checkBtf(A, 4, new long[] {2, 3, 4, 1, 0, 5}, new long[] {2, 0, 3, 1}, new long[] {0, 4, 5, 6, 6},
        new long[] {0, 1, 2, 3, 4}, new long[] {0, 3, 4, 6, 6}, new long[] {0, 1, 3, 3, 4});

    A = new DM(6, 4);
    A.set(0, 0, 1);
    A.set(1, 2, 1);
    A.set(2, 2, 1);
    A.set(5, 3, 1);
    A.set(4, 1, 1);
    A.set(3, 0, 1);
    A.sparsity().spy();
    checkBtf(A, 3, new long[] {2, 3, 0, 1, 4, 5}, new long[] {0, 2, 1, 3}, new long[] {0, 4, 5, 6},
        new long[] {0, 2, 3, 4}, new long[] {0, 2, 4, 6, 6}, new long[] {0, 2, 4, 4, 4});

    A = new DM(6, 4);
    A.set(0, 0, 1);
    A.set(1, 2, 1);
    A.set(2, 2, 1);
    A.set(5, 3, 1);
    A.set(4, 1, 1);
    A.set(3, 0, 1);
    A = A.plus(DM.eye(6).get(new Slice(), new Slice(0, 4)));
    A.sparsity().spy();
    checkBtf(A, 1, new long[] {4, 5, 0, 1, 2, 3}, new long[] {0, 1, 2, 3}, new long[] {0, 6}, new long[] {0, 4},
        new long[] {0, 2, 6, 6, 6}, new long[] {0, 4, 4, 4, 4});
  }

  private static Function jacobian_old(Function f, long i, long j) {
    List<String> out = new ArrayList<>();
    out.add("jac:" + f.name_out(j) + ":" + f.name_in(i));
    out.addAll(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), out);
  }

  @Test
  public void test_jacsparsityHierarchical() {
    SX X = SX.sym("X", 100);
    SX P = SX.sym("P", 1000);
    SX optvar = ca.vertcat(List.of(X, P));
    SX p = SX.sym("p");
    Function g = new Function("g", List.of(optvar, p), List.of(X.times(p)), Map.of("verbose", true));
    Function J = jacobian_old(g, 0, 0);
    assertTrue(new DM(J.sparsity_out(0)).get(new Slice(), new Slice(0, X.nnz())).sparsity()
        .is_equal(Sparsity.diag(100)));

    X = SX.sym("X", 100);
    P = SX.sym("P", 1000);
    p = SX.sym("p");
    g = new Function("g", List.of(X, P, p), List.of(ca.vertcat(List.of(X.times(p), P))), Map.of("verbose", true));
    J = jacobian_old(g, 0, 0);
    assertTrue(new DM(J.sparsity_out(0)).get(new Slice(0, X.nnz()), new Slice()).sparsity()
        .is_equal(Sparsity.diag(100)));
  }

  @Test
  public void test_jacsparsityHierarchicalSymm() {
    GlobalOptions.setHierarchicalSparsity(false);
    try {
      Sparsity sp = Sparsity.banded(4129, 1);
      MX x = MX.sym("x", sp.size1());
      MX H = ca.bilin(MX.ones(sp), x, x);
      Sparsity sp2 = ((MX) ca.hessian(H, x).get(0)).sparsity();
      assertTrue(sp.is_equal(sp2));
    } finally {
      GlobalOptions.setHierarchicalSparsity(true);
    }
  }

  @Test
  public void test_rowcol() {
    long n = 3;
    Sparsity s = Sparsity.rowcol(new long[] {n - 1, 0}, new long[] {0, n - 1}, n, n);
    assertArrayEquals(new long[] {0, 2, 2, 4}, s.colind());
    assertArrayEquals(new long[] {0, 2, 0, 2}, s.row());
  }

  @Test
  public void test_inverse() {
    DM d0 = randDMSymm(20, 20, 0.6);
    for (Sparsity sp : List.of(d0.sparsity(), Sparsity.dense(4, 4), new Sparsity(4, 4), Sparsity.lower(4),
        Sparsity.lower(4).T())) {
      DM d = DM.ones(sp);
      DM dt = ca.sparsify(ca.minus(1, d));
      dt = DM.ones(dt.sparsity());
      DM trial = DM.ones(sp.pattern_inverse());
      d.print_dense();
      dt.print_dense();
      trial.print_dense();
      checkarray(trial, dt);
    }
  }

  @Test
  public void test_kron() {
    DM a = ca.sparsify(dm(new double[][] {{1, 0, 6}, {2, 7, 0}}));
    DM b = ca.sparsify(dm(new double[][] {{1, 0, 0}, {2, 3, 7}, {0, 0, 9}, {1, 12, 13}}));
    Sparsity c_ = ca.kron(a.sparsity(), b.sparsity());
    assertEquals(a.size1() * b.size1(), c_.size1());
    assertEquals(a.size2() * b.size2(), c_.size2());
    assertEquals(a.nnz() * b.nnz(), c_.nnz());
    checkarray(new DM(c_, 1), new DM(ca.kron(a, b).sparsity(), 1));
  }

  @Test
  public void test_nz_method() {
    int n = 20, m = 25;
    DM d = randDM(n, m, 0.5);
    DM D = ca.densify(ca.vec(d));
    DM dn = new DM(d.nonzeros());
    Sparsity sp = d.sparsity();
    TreeSet<Long> zs = new TreeSet<>();
    for (int i = 0; i < 200; i++) zs.add((long) random.nextInt(n * m));
    long[] z = toLongs(new ArrayList<>(zs));
    long[] zres = sp.get_nz(z);
    List<Long> ia = new ArrayList<>(), ib = new ArrayList<>(), ic = new ArrayList<>();
    for (int k = 0; k < z.length; k++) {
      if (zres[k] >= 0) {
        ia.add(zres[k]);
        ib.add(z[k]);
      } else {
        ic.add(z[k]);
      }
    }
    checkarray(dn.get(idx(toLongs(ia))), D.get(idx(toLongs(ib))));
    for (double e : D.get(idx(toLongs(ic))).nonzeros()) assertEquals(0.0, e);
  }

  @Test
  public void test_serialize() {
    for (Sparsity a : List.of(new Sparsity(), Sparsity.dense(4, 5), Sparsity.lower(5))) {
      Sparsity b = Sparsity.deserialize(a.serialize());
      if (a.is_null()) {
        assertTrue(b.is_null());
      } else {
        checkarray(new DM(a, 1), new DM(b, 1));
      }
    }
  }

  @Test
  public void test_is_subset() {
    Sparsity[][] pairs = {{Sparsity.lower(3), Sparsity.dense(3, 3)}, {Sparsity.diag(3), Sparsity.dense(3, 3)},
        {Sparsity.diag(3), Sparsity.lower(3)}, {new Sparsity(3, 3), Sparsity.lower(3)}};
    for (Sparsity[] p : pairs) {
      assertTrue(p[0].is_subset(p[1]));
      assertFalse(p[1].is_subset(p[0]));
    }
  }

  @Test
  public void test_is_compactible() {
    message("is_compactible");
    Sparsity sp = new Sparsity(4, 5);
    sp.add_nz(0, 1);
    sp.add_nz(2, 1);
    sp.add_nz(0, 3);
    sp.add_nz(2, 3);
    List<Object> r = sp.is_compactible();
    assertTrue((Boolean) r.get(0));
    assertArrayEquals(new long[] {0, 2}, longs(r.get(1)));
    assertArrayEquals(new long[] {1, 3}, longs(r.get(2)));

    Sparsity sp_partial = new Sparsity(4, 5);
    sp_partial.add_nz(0, 1);
    sp_partial.add_nz(2, 1);
    sp_partial.add_nz(2, 3);
    assertFalse((Boolean) sp_partial.is_compactible().get(0));

    r = Sparsity.rowcol(new long[] {0, 2}, new long[] {1, 3}, 4, 5).is_compactible();
    assertTrue((Boolean) r.get(0));
    assertArrayEquals(new long[] {0, 2}, longs(r.get(1)));
    assertArrayEquals(new long[] {1, 3}, longs(r.get(2)));

    r = new Sparsity(3, 4).is_compactible();
    assertTrue((Boolean) r.get(0));
    assertArrayEquals(new long[] {}, longs(r.get(1)));
    assertArrayEquals(new long[] {}, longs(r.get(2)));

    r = Sparsity.dense(3, 4).is_compactible();
    assertTrue((Boolean) r.get(0));
    assertArrayEquals(new long[] {0, 1, 2}, longs(r.get(1)));
    assertArrayEquals(new long[] {0, 1, 2, 3}, longs(r.get(2)));

    assertFalse((Boolean) Sparsity.diag(4).is_compactible().get(0));
  }

  @Test
  public void test_star_coloring_cutoff() {
    Sparsity sp = Sparsity.dense(4, 4);
    for (long ordering : new long[] {0, 1}) {
      assertTrue(sp.star_coloring(ordering, 2).is_null());
      assertTrue(sp.star_coloring2(ordering, 2).is_null());
      assertEquals(4, sp.star_coloring(ordering, 4).size2());
    }
    for (boolean new_algo : new boolean[] {false, true}) {
      List<Object> r = sp.star_coloring_new(Map.of("cutoff", 2, "new_algo", new_algo));
      assertTrue(((Sparsity) r.get(0)).is_null());
    }
  }

  @Test
  public void test_compressed_dense_short_form() {
    for (Sparsity sp : List.of(Sparsity.dense(2, 3), Sparsity.dense(2, 0), Sparsity.dense(1, 1))) {
      long[] v = sp.compress(false);
      assertEquals(3, v.length);
      assertTrue(Sparsity.compressed(v).is_equal(sp));
    }
  }

  @Test
  public void test_repr_el_range() {
    Sparsity sp = Sparsity.triplet(3, 4, new long[] {0, 2, 1, 0}, new long[] {0, 0, 2, 3});
    assertEquals("nonzero index 3 (row 0, col 3)", sp.repr_el(3));
    for (long k : new long[] {-1, 4, 5}) assertInException("out of range for 4 nonzeros", () -> sp.repr_el(k));
  }

  @Test
  public void test_get_nz_range() {
    for (Sparsity sp : List.of(Sparsity.dense(2, 2), Sparsity.triplet(3, 4, new long[] {0, 2, 1, 0},
        new long[] {0, 0, 2, 3}))) {
      assertInException("out of bounds", () -> sp.get_nz(new long[] {0, 99}));
    }
    assertArrayEquals(new long[] {3, -1}, Sparsity.dense(2, 2).get_nz(new long[] {3, -1}));
  }

  @Test
  public void test_star_coloring_symmetric() {
    for (Sparsity sp : List.of(Sparsity.upper(3), Sparsity.dense(2, 3))) {
      assertInException("requires a symmetric pattern", () -> sp.star_coloring());
      assertInException("requires a symmetric pattern", () -> sp.star_coloring2());
    }
  }

  @Test
  public void test_banded_wide() {
    assertTrue(Sparsity.banded(4, 5).is_equal(Sparsity.dense(4, 4)));
    assertTrue(Sparsity.banded(4, 3).is_equal(Sparsity.dense(4, 4)));
    assertArrayEquals(new long[] {0, 0}, Sparsity.banded(0, 2).size());
    assertInException("should be smaller than size", () -> Sparsity.band(4, 4));
  }

  @Test
  public void test_density_empty() {
    for (Sparsity sp : List.of(new Sparsity(0, 3), new Sparsity(0, 0), Sparsity.dense(3, 0)))
      assertEquals(100.0, sp.density());
  }

  @Test
  public void test_triplet_bounds_message() {
    assertInException("Row index (2) out of bounds [0,2[",
        () -> Sparsity.triplet(2, 2, new long[] {2}, new long[] {0}));
  }

  @Test
  public void test_scc_square() {
    assertInException("Expected a square pattern", () -> Sparsity.dense(2, 3).scc());
  }

  @Test
  public void test_eq_foreign_operand() {
    // Python data-model parts (== against int/str, TypeError on sp + 3) have no Java counterpart
    Sparsity sp = Sparsity.diag(3);
    assertTrue(sp.is_equal(Sparsity.diag(3)));
  }
}
