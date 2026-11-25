using System;
using System.Collections.Generic;
using System.Linq;

namespace CasADi.Tests {

/// <summary>Port of test/java SparsityTest (test/python/sparsity.py)</summary>
public class SparsityTest : CasadiTestCase {
  static Sparsity fromNz(long n, long m, long[,] nz) {
    var a = ca.Sparsity(n, m);
    for (int k = 0; k < nz.GetLength(0); k++) a.add_nz(nz[k, 0], nz[k, 1]);
    return a;
  }

  static HashSet<(long, long)> nzSet(long[,] nz) {
    var r = new HashSet<(long, long)>();
    for (int k = 0; k < nz.GetLength(0); k++) r.Add((nz[k, 0], nz[k, 1]));
    return r;
  }

  static DM idx(IEnumerable<long> i) => i.ToArray();

  static void assertSeq(IEnumerable<long> expected, IEnumerable<long> actual) =>
    Check.True(expected.SequenceEqual(actual),
               "expected [" + string.Join(", ", expected) + "], got [" + string.Join(", ", actual) + "]");

  /// <summary>Python helper randDM(n, m, sparsity, symm=True)</summary>
  DM randDMSymm(long n, long m, double density) {
    DM d = density < 1 ? randDM(n, m, density) : randDM(n, m);
    return (d + d.T()) / 2.0;
  }

  public void test_union() {
    message("Sparsity union");
    long[,] nza = { { 0, 0 }, { 0, 1 }, { 2, 0 }, { 3, 1 } };
    long[,] nzb = { { 0, 2 }, { 0, 0 }, { 2, 2 } };
    Sparsity a = fromNz(4, 5, nza), b = fromNz(4, 5, nzb);
    Sparsity c = a.unite(b);
    c = a.plus(b);
    var u = nzSet(nza);
    u.UnionWith(nzSet(nzb));
    assertEquals((long)u.Count, c.nnz());
    var col = c.get_col();
    for (int k = 0; k < c.nnz(); k++) assertTrue(u.Contains((c.row(k), col[k])));
  }

  public void test_intersection() {
    message("Sparsity intersection");
    long[,] nza = { { 0, 0 }, { 0, 1 }, { 2, 0 }, { 3, 1 }, { 2, 3 } };
    long[,] nzb = { { 0, 2 }, { 0, 0 }, { 2, 2 }, { 2, 3 } };
    Sparsity a = fromNz(4, 5, nza), b = fromNz(4, 5, nzb);
    var sa = nzSet(nza);
    var sb = nzSet(nzb);
    Sparsity c = a.intersect(b);
    for (int k = 0; k < c.nnz(); k++) {
      var ind = (c.row(k), c.get_col()[k]);
      assertTrue(sa.Contains(ind) && sb.Contains(ind));
    }
    c = a.times(b);
    var i = new HashSet<(long, long)>(sa);
    i.IntersectWith(sb);
    assertEquals((long)i.Count, c.nnz());
    for (int k = 0; k < c.nnz(); k++) {
      var ind = (c.row(k), c.get_col()[k]);
      assertTrue(sa.Contains(ind) && sb.Contains(ind));
    }
  }

  public void test_get_nzDense() {
    message("get_nzDense");
    Sparsity a = fromNz(4, 5, new long[,] { { 0, 0 }, { 0, 1 }, { 2, 0 }, { 3, 1 } });
    DM A = DM.ones(a);
    DM Ad = A.ToArray();
    foreach (long i in a.find()) assertEquals(1.0, (double)Ad.get_nz(false, ca.DM((double)i)));
  }

  public void test_find_nonzero() {
    DM d = randDM(20, 10, 0.6);
    Sparsity sp = d.sparsity();
    Sparsity sp2 = Sparsity.nonzeros(20, 10, sp.find());
    assertTrue(sp.is_equal(sp2));
  }

  public void test_enlarge() {
    message("enlarge");
    message(":dense");
    Sparsity sp = Sparsity.dense(3, 4);
    long[] col = { 1, 2, 4 };
    long[] row = { 0, 3, 4, 6 };
    sp.enlarge(7, 8, col, row);
    var z = new double[7, 8];
    foreach (long i in col)
      foreach (long j in row) z[i, j] = 1;
    checkarray(DM.ones(sp), z);
    message(":sparse");
    sp = ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 }).T();
    DM n = DM.ones(sp);
    z = new double[7, 8];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++) z[col[i], row[j]] = (double)n[i, j];
    sp.enlarge(7, 8, new long[] { 1, 2, 4 }, new long[] { 0, 3, 4, 6 });
    checkarray(DM.ones(sp), z);
  }

  static DM tomatrix(Sparsity s) {
    DM d = DM.ones(s);
    for (long k = 0; k < d.nnz(); k++) d.set_nz(k + 1, false, ca.DM((double)k));
    return d;
  }

  public void test_NZ() {
    message("NZ constructor");
    long[,] nza = { { 0, 0 }, { 0, 1 }, { 2, 0 }, { 2, 3 }, { 2, 4 }, { 3, 1 } };
    Sparsity a = fromNz(4, 5, nza);
    var r = new long[nza.GetLength(0)];
    var c = new long[nza.GetLength(0)];
    for (int k = 0; k < r.Length; k++) {
      r[k] = nza[k, 0];
      c[k] = nza[k, 1];
    }
    Sparsity b = Sparsity.triplet(4, 5, r, c);
    checkarray(tomatrix(a), tomatrix(b), "rowcol");
  }

  /// <summary>First test_rowcol of sparsity.py, shadowed there by the later definition</summary>
  public void test_rowcol_shadowed() {
    message("rowcol constructor");
    long[] r = { 0, 1, 3 };
    long[] c = { 1, 4 };
    var a = ca.Sparsity(4, 5);
    foreach (long i in r)
      foreach (long j in c) a.add_nz(i, j);
    Sparsity b = Sparsity.rowcol(r, c, 4, 5);
    checkarray(tomatrix(a), tomatrix(b), "rowcol");
  }

  public void test_reshape() {
    message("Reshape");
    Sparsity a = fromNz(4, 5, new long[,] { { 0, 0 }, { 0, 1 }, { 2, 0 }, { 2, 3 }, { 2, 4 }, { 3, 1 } });
    double[,] A = tomatrix(a).ToArray();
    DM B = tomatrix(ca.reshape(a, 2, 10));
    var B_ = new double[2, 10];
    for (int k = 0; k < 20; k++) B_[k % 2, k / 2] = A[k % 4, k / 4];
    checkarray(B, B_);
  }

  // test_vec: returns immediately in Python ("doesn't make much sense")

  public void test_refcount() {
    var x = ca.DM(Sparsity.lower(4), 5);
    Sparsity s = ca.mtimes(x, x).sparsity();
    assertEquals(16L, s.numel());
  }

  public void test_splower() {
    var sp = ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 });
    Console.WriteLine(ca.DM(sp, 1));
    Console.WriteLine(ca.DM(ca.tril(sp), 1));
    Console.WriteLine("[" + string.Join(", ", sp.get_lower()) + "]");
  }

  public void test_diag() {
    message("diag");
    var A = ca.Sparsity(5, 5);
    A.add_nz(1, 1);
    A.add_nz(2, 4);
    A.add_nz(3, 3);
    var (s, m) = A.get_diag();
    checkarray(DM.ones(s), new double[,] { { 0 }, { 1 }, { 0 }, { 1 }, { 0 } });
    assertSeq(new long[] { 0, 1 }, m);

    A = ca.Sparsity(5, 1);
    A.add_nz(1, 0);
    A.add_nz(2, 0);
    A.add_nz(4, 0);
    double[,] reference = { { 0, 0, 0, 0, 0 }, { 0, 1, 0, 0, 0 }, { 0, 0, 1, 0, 0 }, { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 1 } };
    (s, m) = A.get_diag();
    checkarray(DM.ones(s), reference);
    assertSeq(new long[] { 0, 1, 2 }, m);

    A = ca.Sparsity(1, 5);
    A.add_nz(0, 1);
    A.add_nz(0, 2);
    A.add_nz(0, 4);
    (s, m) = A.get_diag();
    checkarray(DM.ones(s), reference);
    assertSeq(new long[] { 0, 1, 2 }, m);
  }

  public void test_sparsityindex() {
    message("sparsity indexing");
    long[,] nza = { { 0, 0 }, { 0, 1 }, { 2, 0 }, { 2, 3 }, { 3, 3 }, { 2, 4 }, { 3, 1 }, { 4, 1 } };
    Sparsity a = fromNz(5, 5, nza);
    MX b = MX.sym("b", a);
    Check.Throws<Exception>(() => b.get(Sparsity.diag(3)));
    Sparsity d = Sparsity.diag(5);
    MX c = b.get(d);
    assertTrue(c.sparsity().is_equal(d));
    var f = ca.Function("f", new[] { b }, new[] { c });
    var v = new double[nza.GetLength(0)];
    for (int k = 0; k < v.Length; k++) v[k] = k + 1;
    var fin = ca.DM(b.sparsity(), v);
    DM f_out = f.Call(fin);
    checkarray(f_out.nonzeros().ToArray(), ca.vertcat(1, 0, 0, 7, 0), "sparsity index");
  }

  public void test_get_ccs() {
    message("CCS format");
    Sparsity a = fromNz(4, 5, new long[,] { { 0, 0 }, { 0, 1 }, { 2, 0 }, { 2, 3 }, { 3, 3 }, { 2, 4 }, { 3, 1 } });
    var (c1, r1) = a.get_ccs();
    var (c2, r2) = a.T().get_crs();
    Console.WriteLine("[" + string.Join(", ", c1) + "] [" + string.Join(", ", r1) + "]");
    Console.WriteLine("[" + string.Join(", ", c2) + "] [" + string.Join(", ", r2) + "]");
  }

  List<long> randomPermutation(long n) {
    var p = new List<long>();
    for (long i = 0; i < n; i++) p.Add(i);
    for (int i = p.Count - 1; i > 0; i--) {
      int j = random.Next(i + 1);
      (p[i], p[j]) = (p[j], p[i]);
    }
    return p;
  }

  static List<long> blocks(DM Ar) {
    Sparsity ST = Ar.sparsity();
    var blocks = new List<long>();
    long acc = -1, mc = 0;
    List<long> row = ST.row(), colind = ST.colind();
    for (int i = 0; i < Ar.size1(); i++) {
      mc = Math.Max(row[(int)colind[i + 1] - 1], mc);
      if (mc == i) {
        blocks.Add(i - acc);
        acc = i;
      }
    }
    return blocks;
  }

  public void test_dm_diagcat_dense() {
    message("Dulmage-Mendelsohn");
    for (int k = 0; k < 20; k++) {
      var Ai = new List<DM>();
      for (int j = 0; j < 10; j++) {
        long d = 1 + random.Next(10);
        Ai.Add(randDM(d, d));
      }
      DM A = ca.diagcat(Ai);
      DM perm = idx(randomPermutation(A.size1()));
      DM AP = A.get(perm, perm);
      var r = AP.sparsity().btf();
      DM Ar = AP.get(idx(r.Item2), idx(r.Item3));
      var truth = Ai.Select(a => a.size1()).OrderBy(x => x).ToList();
      var tryme = blocks(Ar).OrderBy(x => x).ToList();
      assertSeq(truth, tryme);
    }
  }

  public void test_scc_diagcat_sparse() {
    message("scc");
    for (int k = 0; k < 20; k++) {
      var Ai = new List<DM>();
      for (int j = 0; j < 10; j++) {
        long d = 1 + random.Next(10);
        Ai.Add(randDMSymm(d, d, 0.6));
      }
      DM A = ca.diagcat(Ai);
      DM perm = idx(randomPermutation(A.size1()));
      DM AP = A.get(perm, perm);
      var (n, pv, _) = AP.sparsity().scc();
      DM p = idx(pv);
      DM Ar = AP.get(p, p);
      var tryme = blocks(Ar);
      assertTrue(n >= Ai.Count);
      assertTrue(n >= tryme.Count);
    }
  }

  static void checkBtf(DM A, long ret, long[] rowperm, long[] colperm, long[] rowblock, long[] colblock,
                       long[] coarse_rowblock, long[] coarse_colblock) {
    var r = A.sparsity().btf();
    assertEquals(ret, r.Item1);
    assertSeq(rowperm, r.Item2);
    assertSeq(colperm, r.Item3);
    assertSeq(rowblock, r.Item4);
    assertSeq(colblock, r.Item5);
    assertSeq(coarse_rowblock, r.Item6);
    assertSeq(coarse_colblock, r.Item7);
  }

  public void test_dm() {
    var A = ca.DM(6, 4);
    A[0, 0] = 1;
    A[1, 2] = 1;
    A[2, 2] = 1;
    A[5, 3] = 1;
    checkBtf(A, 4, new long[] { 2, 3, 4, 1, 0, 5 }, new long[] { 2, 0, 3, 1 }, new long[] { 0, 4, 5, 6, 6 },
             new long[] { 0, 1, 2, 3, 4 }, new long[] { 0, 3, 4, 6, 6 }, new long[] { 0, 1, 3, 3, 4 });

    A = ca.DM(6, 4);
    A[0, 0] = 1;
    A[1, 2] = 1;
    A[2, 2] = 1;
    A[5, 3] = 1;
    A[4, 1] = 1;
    A[3, 0] = 1;
    A.sparsity().spy();
    checkBtf(A, 3, new long[] { 2, 3, 0, 1, 4, 5 }, new long[] { 0, 2, 1, 3 }, new long[] { 0, 4, 5, 6 },
             new long[] { 0, 2, 3, 4 }, new long[] { 0, 2, 4, 6, 6 }, new long[] { 0, 2, 4, 4, 4 });

    A = ca.DM(6, 4);
    A[0, 0] = 1;
    A[1, 2] = 1;
    A[2, 2] = 1;
    A[5, 3] = 1;
    A[4, 1] = 1;
    A[3, 0] = 1;
    A = A + DM.eye(6)[.., 0..4];
    A.sparsity().spy();
    checkBtf(A, 1, new long[] { 4, 5, 0, 1, 2, 3 }, new long[] { 0, 1, 2, 3 }, new long[] { 0, 6 }, new long[] { 0, 4 },
             new long[] { 0, 2, 6, 6, 6 }, new long[] { 0, 4, 4, 4, 4 });
  }

  static Function jacobian_old(Function f, long i, long j) {
    var outs = new List<string> { "jac:" + f.name_out(j) + ":" + f.name_in(i) };
    outs.AddRange(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), outs);
  }

  public void test_jacsparsityHierarchical() {
    SX X = SX.sym("X", 100);
    SX P = SX.sym("P", 1000);
    SX optvar = ca.vertcat(X, P);
    SX p = SX.sym("p");
    var verbose = new Dictionary<string, object> { ["verbose"] = true };
    var g = ca.Function("g", new[] { optvar, p }, new[] { X * p }, verbose);
    Function J = jacobian_old(g, 0, 0);
    assertTrue(ca.DM(J.sparsity_out(0))[.., 0..(int)X.nnz()].sparsity().is_equal(Sparsity.diag(100)));

    X = SX.sym("X", 100);
    P = SX.sym("P", 1000);
    p = SX.sym("p");
    g = ca.Function("g", new[] { X, P, p }, new[] { ca.vertcat(X * p, P) }, verbose);
    J = jacobian_old(g, 0, 0);
    assertTrue(ca.DM(J.sparsity_out(0))[0..(int)X.nnz(), ..].sparsity().is_equal(Sparsity.diag(100)));
  }

  public void test_jacsparsityHierarchicalSymm() {
    GlobalOptions.setHierarchicalSparsity(false);
    try {
      Sparsity sp = Sparsity.banded(4129, 1);
      MX x = MX.sym("x", sp.size1());
      MX H = ca.bilin(MX.ones(sp), x, x);
      Sparsity sp2 = ca.hessian(H, x).Item1.sparsity();
      assertTrue(sp.is_equal(sp2));
    } finally {
      GlobalOptions.setHierarchicalSparsity(true);
    }
  }

  public void test_rowcol() {
    long n = 3;
    Sparsity s = Sparsity.rowcol(new long[] { n - 1, 0 }, new long[] { 0, n - 1 }, n, n);
    assertSeq(new long[] { 0, 2, 2, 4 }, s.colind());
    assertSeq(new long[] { 0, 2, 0, 2 }, s.row());
  }

  public void test_inverse() {
    DM d0 = randDMSymm(20, 20, 0.6);
    foreach (Sparsity sp in new[] { d0.sparsity(), Sparsity.dense(4, 4), ca.Sparsity(4, 4), Sparsity.lower(4),
                                    Sparsity.lower(4).T() }) {
      DM d = DM.ones(sp);
      DM dt = ca.sparsify(1 - d);
      dt = DM.ones(dt.sparsity());
      DM trial = DM.ones(sp.pattern_inverse());
      d.print_dense();
      dt.print_dense();
      trial.print_dense();
      checkarray(trial, dt);
    }
  }

  public void test_kron() {
    DM a = ca.sparsify(ca.DM(new double[,] { { 1, 0, 6 }, { 2, 7, 0 } }));
    DM b = ca.sparsify(ca.DM(new double[,] { { 1, 0, 0 }, { 2, 3, 7 }, { 0, 0, 9 }, { 1, 12, 13 } }));
    Sparsity c_ = ca.kron(a.sparsity(), b.sparsity());
    assertEquals(a.size1() * b.size1(), c_.size1());
    assertEquals(a.size2() * b.size2(), c_.size2());
    assertEquals(a.nnz() * b.nnz(), c_.nnz());
    checkarray(ca.DM(c_, 1), ca.DM(ca.kron(a, b).sparsity(), 1));
  }

  public void test_nz_method() {
    int n = 20, m = 25;
    DM d = randDM(n, m, 0.5);
    DM D = ca.densify(ca.vec(d));
    DM dn = d.nonzeros().ToArray();
    Sparsity sp = d.sparsity();
    var zs = new SortedSet<long>();
    for (int i = 0; i < 200; i++) zs.Add(random.Next(n * m));
    long[] z = zs.ToArray();
    var zres = sp.get_nz(z);
    List<long> ia = new List<long>(), ib = new List<long>(), ic = new List<long>();
    for (int k = 0; k < z.Length; k++) {
      if (zres[k] >= 0) {
        ia.Add(zres[k]);
        ib.Add(z[k]);
      } else {
        ic.Add(z[k]);
      }
    }
    checkarray(dn[ia.ToArray()], D[ib.ToArray()]);
    foreach (double e in D[ic.ToArray()].nonzeros()) assertEquals(0.0, e);
  }

  public void test_serialize() {
    foreach (Sparsity a in new[] { ca.Sparsity(), Sparsity.dense(4, 5), Sparsity.lower(5) }) {
      Sparsity b = Sparsity.deserialize(a.serialize());
      if (a.is_null()) {
        assertTrue(b.is_null());
      } else {
        checkarray(ca.DM(a, 1), ca.DM(b, 1));
      }
    }
  }

  public void test_is_subset() {
    var pairs = new[] { (Sparsity.lower(3), Sparsity.dense(3, 3)), (Sparsity.diag(3), Sparsity.dense(3, 3)),
                        (Sparsity.diag(3), Sparsity.lower(3)), (ca.Sparsity(3, 3), Sparsity.lower(3)) };
    foreach (var (p0, p1) in pairs) {
      assertTrue(p0.is_subset(p1));
      assertFalse(p1.is_subset(p0));
    }
  }

  public void test_is_compactible() {
    message("is_compactible");
    var sp = ca.Sparsity(4, 5);
    sp.add_nz(0, 1);
    sp.add_nz(2, 1);
    sp.add_nz(0, 3);
    sp.add_nz(2, 3);
    var (ok, rr, cc) = sp.is_compactible();
    assertTrue(ok);
    assertSeq(new long[] { 0, 2 }, rr);
    assertSeq(new long[] { 1, 3 }, cc);

    var sp_partial = ca.Sparsity(4, 5);
    sp_partial.add_nz(0, 1);
    sp_partial.add_nz(2, 1);
    sp_partial.add_nz(2, 3);
    assertFalse(sp_partial.is_compactible().Item1);

    (ok, rr, cc) = Sparsity.rowcol(new long[] { 0, 2 }, new long[] { 1, 3 }, 4, 5).is_compactible();
    assertTrue(ok);
    assertSeq(new long[] { 0, 2 }, rr);
    assertSeq(new long[] { 1, 3 }, cc);

    (ok, rr, cc) = ca.Sparsity(3, 4).is_compactible();
    assertTrue(ok);
    assertSeq(new long[] { }, rr);
    assertSeq(new long[] { }, cc);

    (ok, rr, cc) = Sparsity.dense(3, 4).is_compactible();
    assertTrue(ok);
    assertSeq(new long[] { 0, 1, 2 }, rr);
    assertSeq(new long[] { 0, 1, 2, 3 }, cc);

    assertFalse(Sparsity.diag(4).is_compactible().Item1);
  }

  public void test_star_coloring_cutoff() {
    Sparsity sp = Sparsity.dense(4, 4);
    foreach (long ordering in new long[] { 0, 1 }) {
      assertTrue(sp.star_coloring(ordering, 2).is_null());
      assertTrue(sp.star_coloring2(ordering, 2).is_null());
      assertEquals(4L, sp.star_coloring(ordering, 4).size2());
    }
    foreach (bool new_algo in new[] { false, true }) {
      var r = sp.star_coloring_new(new Dictionary<string, object> { ["cutoff"] = 2L, ["new_algo"] = new_algo });
      assertTrue(r.Item1.is_null());
    }
  }

  public void test_compressed_dense_short_form() {
    foreach (Sparsity sp in new[] { Sparsity.dense(2, 3), Sparsity.dense(2, 0), Sparsity.dense(1, 1) }) {
      var v = sp.compress(false);
      assertEquals(3, v.Count);
      assertTrue(Sparsity.compressed(v).is_equal(sp));
    }
  }

  public void test_repr_el_range() {
    Sparsity sp = Sparsity.triplet(3, 4, new long[] { 0, 2, 1, 0 }, new long[] { 0, 0, 2, 3 });
    assertEquals("nonzero index 3 (row 0, col 3)", sp.repr_el(3));
    foreach (long k in new long[] { -1, 4, 5 }) assertInException("out of range for 4 nonzeros", () => sp.repr_el(k));
  }

  public void test_get_nz_range() {
    foreach (Sparsity sp in new[] { Sparsity.dense(2, 2),
                                    Sparsity.triplet(3, 4, new long[] { 0, 2, 1, 0 }, new long[] { 0, 0, 2, 3 }) }) {
      assertInException("out of bounds", () => sp.get_nz(new long[] { 0, 99 }));
    }
    assertSeq(new long[] { 3, -1 }, Sparsity.dense(2, 2).get_nz(new long[] { 3, -1 }));
  }

  public void test_star_coloring_symmetric() {
    foreach (Sparsity sp in new[] { Sparsity.upper(3), Sparsity.dense(2, 3) }) {
      assertInException("requires a symmetric pattern", () => sp.star_coloring());
      assertInException("requires a symmetric pattern", () => sp.star_coloring2());
    }
  }

  public void test_banded_wide() {
    assertTrue(Sparsity.banded(4, 5).is_equal(Sparsity.dense(4, 4)));
    assertTrue(Sparsity.banded(4, 3).is_equal(Sparsity.dense(4, 4)));
    assertEquals((0L, 0L), Sparsity.banded(0, 2).size());
    assertInException("should be smaller than size", () => Sparsity.band(4, 4));
  }

  public void test_density_empty() {
    foreach (Sparsity sp in new[] { ca.Sparsity(0, 3), ca.Sparsity(0, 0), Sparsity.dense(3, 0) })
      assertEquals(100.0, sp.density());
  }

  public void test_triplet_bounds_message() {
    assertInException("Row index (2) out of bounds [0,2[",
                      () => Sparsity.triplet(2, 2, new long[] { 2 }, new long[] { 0 }));
  }

  public void test_scc_square() {
    assertInException("Expected a square pattern", () => Sparsity.dense(2, 3).scc());
  }

  public void test_eq_foreign_operand() {
    // Python data-model parts (== against int/str, TypeError on sp + 3) have no C# counterpart
    Sparsity sp = Sparsity.diag(3);
    assertTrue(sp.is_equal(Sparsity.diag(3)));
  }
}

}
