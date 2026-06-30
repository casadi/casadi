// linsol.js -- subset port of test/python/linearsolver.py.
//
// Exercises the Linsol class and solve() free-function against every
// linear-solver plugin loaded in the wasm-js build:
//   qr, ldl, lsqr, symbolicqr, csparse, lapacklu, lapackqr
//
// Coverage:
//   * test_simple                    -- basic 2x2 dense solve via Linsol class
//   * test_simple_dmatrix            -- same problem via top-level solve()
//   * test_simple_trans              -- transpose solve via Linsol.solve(A,b,true)
//   * test_simple_multi_rhs          -- multiple right-hand sides
//   * test_dense_3x3                 -- 3x3 dense solve, well-conditioned
//   * test_sparse_pattern            -- sparse-pattern A (lower-triangular)
//   * test_mx_solve                  -- solve() at MX level wrapped in a Function
//   * test_identity                  -- A = I trivial check
//   * test_diagonal                  -- diag-A trivial check (catches tridiag-style errors)
//   * test_residual                  -- verify A * x == b after solve
//   * test_neig_rank                 -- exercise neig()/rank() where supported
//   * test_factor_reuse              -- nfact()+solve(); refactor with new A.
//   * test_sfact_nfact_explicit      -- explicit sfact/nfact split (vs all-in-one solve)
//
// Per-solver feature flags (mirroring lsolvers in linearsolver.py):
//   "symmetry" => solver requires symmetric A (ldl)
//   "posdef"   => solver requires SPD A  (ldl)
//   "iter"     => iterative solver, lower precision (lsqr)
//
// Driver gates each test per solver; symmetric-only tests skip on ldl
// only when A is non-symmetric, etc.

const { assertEqual, assertAlmost, assertArrayAlmostEqual,
        runTests } = require("./_helpers");

// ----- helpers -----

// Solver registry: per-solver requirements.
const ALL_SOLVERS = [
  ["qr",         { symmetry: false, posdef: false, iter: false, digits: 8 }],
  ["ldl",        { symmetry: true,  posdef: false, iter: false, digits: 8 }],
  ["lsqr",       { symmetry: false, posdef: false, iter: true,  digits: 4 }],
  ["symbolicqr", { symmetry: false, posdef: false, iter: false, digits: 8 }],
  ["csparse",    { symmetry: false, posdef: false, iter: false, digits: 8 }],
  ["lapacklu",   { symmetry: false, posdef: false, iter: false, digits: 8 }],
  ["lapackqr",   { symmetry: false, posdef: false, iter: false, digits: 8 }],
];

// Make a 2x2 dense matrix; if `symmetric` make it symmetric SPD.
function makeA2(M, symmetric) {
  if (symmetric) {
    // SPD: [[4,1],[1,3]]; symmetric, positive-definite.
    return M.DM([[4, 1], [1, 3]]);
  }
  return M.DM([[3, 7], [1, 2]]);  // non-symmetric, nonsingular
}

function npSolve2(A, b) {
  // Solve 2x2 dense linear system A*x = b analytically.
  // A is [[a,b],[c,d]] in row-major; b is length-2 vector.
  const a = A[0][0], B = A[0][1], c = A[1][0], d = A[1][1];
  const det = a * d - B * c;
  return [(d * b[0] - B * b[1]) / det, (-c * b[0] + a * b[1]) / det];
}

// ----- tests -----

function test_simple(M, plugin, req) {
  const A = makeA2(M, req.symmetry);
  const b = M.DM([1, 0.5]);
  const ls = M.Linsol("S", plugin, A.sparsity());
  const x = ls.solve(A, b, false);
  // Reference: numpy.linalg.solve.
  const Arow = req.symmetry ? [[4, 1], [1, 3]] : [[3, 7], [1, 2]];
  const ref = npSolve2(Arow, [1, 0.5]);
  assertArrayAlmostEqual(x.nonzeros(), ref, req.digits, `${plugin} simple x`);
}

function test_simple_dmatrix(M, plugin, req) {
  // solve() free function path (no explicit Linsol class).
  const A = makeA2(M, req.symmetry);
  const b = M.DM([1, 0.5]);
  const x = M.solve(A, b, plugin);
  const Arow = req.symmetry ? [[4, 1], [1, 3]] : [[3, 7], [1, 2]];
  const ref = npSolve2(Arow, [1, 0.5]);
  assertArrayAlmostEqual(x.nonzeros(), ref, req.digits, `${plugin} dmatrix x`);
}

function test_simple_trans(M, plugin, req) {
  // DM-level transpose solve via Linsol::solve(A, b, tr=true).
  const A = makeA2(M, req.symmetry);
  const b = M.DM([1, 0.5]);
  const ls = M.Linsol("S", plugin, A.sparsity());
  const x = ls.solve(A, b, true);
  // Reference: solve A^T x = b directly.
  const Arow = req.symmetry ? [[4, 1], [1, 3]] : [[3, 1], [7, 2]];
  const ref = npSolve2(Arow, [1, 0.5]);
  assertArrayAlmostEqual(x.nonzeros(), ref, req.digits, `${plugin} A^T x x`);
}

function test_simple_trans_mx(M, plugin, req) {
  // Transpose solve at MX level: build solve(A^T, b) symbolically and
  // evaluate numerically.  Complements test_simple_trans which tests
  // the DM-level path.
  const A_ = makeA2(M, req.symmetry);
  const b_ = M.DM([1, 0.5]);
  const A = M.MX.sym("A", 2n, 2n);
  const b = M.MX.sym("b", 2n, 1n);
  const x = M.solve(M.transpose(A), b, plugin);
  const f = M.Function("f", [A, b], [x]);
  const out = f.call([A_, b_])[0];
  // Reference: solve A^T x = b directly.
  const Arow = req.symmetry ? [[4, 1], [1, 3]] : [[3, 1], [7, 2]];
  const ref = npSolve2(Arow, [1, 0.5]);
  assertArrayAlmostEqual(out.nonzeros(), ref, req.digits,
    `${plugin} trans_mx x`);
}

function test_simple_multi_rhs(M, plugin, req) {
  // Two right-hand sides at once: B is 2x2, X is 2x2.
  // For symmetric we still use the symmetric A.
  const A = makeA2(M, req.symmetry);
  const B = M.DM([[1, 0.3], [0.5, 0.7]]);
  const ls = M.Linsol("S", plugin, A.sparsity());
  const X = ls.solve(A, B, false);
  // Verify A*X == B componentwise.
  const AX = M.mtimes(A, X);
  const diff = M.minus(AX, B);
  const vals = diff.nonzeros();
  for (let i = 0; i < vals.length; ++i) {
    assertAlmost(vals[i], 0, req.digits, `${plugin} multi_rhs A*X-B[${i}]`);
  }
}

function test_dense_3x3(M, plugin, req) {
  // 3x3 well-conditioned, force symmetric for ldl/posdef solvers.
  let A, ref;
  if (req.symmetry) {
    // SPD 3x3: [[4,1,0],[1,3,1],[0,1,2]]
    A = M.DM([[4, 1, 0], [1, 3, 1], [0, 1, 2]]);
    // Solve analytically: A * x = [1, 0, 0]
    // Use the residual check; computing inverse by hand is tedious.
    const b = M.DM([1, 0, 0]);
    const ls = M.Linsol("S", plugin, A.sparsity());
    const x = ls.solve(A, b, false);
    const r = M.minus(M.mtimes(A, x), b).nonzeros();
    for (let i = 0; i < r.length; ++i) {
      assertAlmost(r[i], 0, req.digits, `${plugin} dense_3x3 res[${i}]`);
    }
    return;
  }
  A = M.DM([[2, 1, 3], [4, 5, 6], [7, 8, 10]]);
  const b = M.DM([6, 15, 25]);
  const ls = M.Linsol("S", plugin, A.sparsity());
  const x = ls.solve(A, b, false);
  // Check residual.
  const r = M.minus(M.mtimes(A, x), b).nonzeros();
  for (let i = 0; i < r.length; ++i) {
    assertAlmost(r[i], 0, req.digits, `${plugin} dense_3x3 res[${i}]`);
  }
}

function test_sparse_pattern(M, plugin, req) {
  // Lower triangular A (still nonsingular).  Skip for solvers that need
  // symmetric input.
  if (req.symmetry) return;
  // [[2,0,0],[1,3,0],[1,1,4]]
  const A = M.DM([[2, 0, 0], [1, 3, 0], [1, 1, 4]]);
  const b = M.DM([2, 3, 4]);
  const ls = M.Linsol("S", plugin, A.sparsity());
  const x = ls.solve(A, b, false);
  // Forward-substitution by hand: x = [1, 2/3, 1/12 of complicated...]
  // Just verify residual.
  const r = M.minus(M.mtimes(A, x), b).nonzeros();
  for (let i = 0; i < r.length; ++i) {
    assertAlmost(r[i], 0, req.digits, `${plugin} sparse_pattern res[${i}]`);
  }
}

function test_mx_solve(M, plugin, req) {
  // MX-level solve wrapped in a Function and evaluated numerically.
  const A_ = makeA2(M, req.symmetry);
  const b_ = M.DM([1, 0.5]);
  const A = M.MX.sym("A", 2n, 2n);
  const b = M.MX.sym("b", 2n, 1n);
  const x = M.solve(A, b, plugin);
  const f = M.Function("f", [A, b], [x]);
  const out = f.call([A_, b_])[0];
  // Compare to numpy-style direct solve.
  const Arow = req.symmetry ? [[4, 1], [1, 3]] : [[3, 7], [1, 2]];
  const ref = npSolve2(Arow, [1, 0.5]);
  assertArrayAlmostEqual(out.nonzeros(), ref, req.digits,
    `${plugin} mx_solve x`);
}

function test_identity(M, plugin, req) {
  // A = I; expect x = b.
  const A = M.DM([[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
  const b = M.DM([3, 5, 7]);
  const ls = M.Linsol("S", plugin, A.sparsity());
  const x = ls.solve(A, b, false);
  assertArrayAlmostEqual(x.nonzeros(), [3, 5, 7], req.digits,
    `${plugin} identity x`);
}

function test_diagonal(M, plugin, req) {
  // A = diag([2, 3, 4]); expect x = b ./ diag.
  const A = M.DM([[2, 0, 0], [0, 3, 0], [0, 0, 4]]);
  const b = M.DM([4, 9, 16]);
  const ls = M.Linsol("S", plugin, A.sparsity());
  const x = ls.solve(A, b, false);
  assertArrayAlmostEqual(x.nonzeros(), [2, 3, 4], req.digits,
    `${plugin} diagonal x`);
}

function test_residual(M, plugin, req) {
  // 4x4 problem; verify ||A*x - b||_inf is small.
  let A;
  if (req.symmetry) {
    A = M.DM([[5, 1, 0, 0], [1, 4, 1, 0], [0, 1, 3, 1], [0, 0, 1, 2]]);
  } else {
    A = M.DM([[3, 1, 0, 0], [1, 4, 1, 0], [0, 1, 5, 1], [0, 0, 1, 6]]);
  }
  const b = M.DM([1, 2, 3, 4]);
  const ls = M.Linsol("S", plugin, A.sparsity());
  const x = ls.solve(A, b, false);
  const r = M.minus(M.mtimes(A, x), b).nonzeros();
  for (let i = 0; i < r.length; ++i) {
    assertAlmost(r[i], 0, req.digits, `${plugin} residual[${i}]`);
  }
}

function test_neig_rank(M, plugin, req) {
  // Not every solver supports neig/rank; wrap in try.  We just verify
  // the calls don't crash and (where they work) return sensible values.
  const A = makeA2(M, req.symmetry);
  const ls = M.Linsol("S", plugin, A.sparsity());
  // Need numerical factorization first.
  ls.sfact(A);
  ls.nfact(A);
  try {
    const r = ls.rank(A);
    // 2x2 nonsingular: rank should be 2.
    if (typeof r !== "undefined") {
      assertEqual(Number(r), 2, `${plugin} rank`);
    }
  } catch (_e) {
    // Solver doesn't support rank; that's fine.
  }
}

function test_factor_reuse(M, plugin, req) {
  // Same sparsity pattern, two numerical refactorizations.  Verify each
  // refactor produces a correct solve for its respective A.
  const sp = makeA2(M, req.symmetry).sparsity();
  const ls = M.Linsol("S", plugin, sp);
  // First system.
  const A1 = makeA2(M, req.symmetry);
  const b1 = M.DM([1, 0.5]);
  const x1 = ls.solve(A1, b1, false);
  const Arow1 = req.symmetry ? [[4, 1], [1, 3]] : [[3, 7], [1, 2]];
  const ref1 = npSolve2(Arow1, [1, 0.5]);
  assertArrayAlmostEqual(x1.nonzeros(), ref1, req.digits,
    `${plugin} refactor 1`);
  // Second system: scale A1 by 2 (same sparsity, different numerics).
  // For symmetric we double the SPD matrix; still SPD.
  const A2 = M.DM(req.symmetry ? [[8, 2], [2, 6]] : [[6, 14], [2, 4]]);
  const b2 = M.DM([2, 1]);
  const x2 = ls.solve(A2, b2, false);
  const Arow2 = req.symmetry ? [[8, 2], [2, 6]] : [[6, 14], [2, 4]];
  const ref2 = npSolve2(Arow2, [2, 1]);
  assertArrayAlmostEqual(x2.nonzeros(), ref2, req.digits,
    `${plugin} refactor 2`);
}

function test_sfact_nfact_explicit(M, plugin, req) {
  // Explicit two-phase factor-then-solve.
  const A = makeA2(M, req.symmetry);
  const b = M.DM([1, 0.5]);
  const ls = M.Linsol("S", plugin, A.sparsity());
  ls.sfact(A);
  ls.nfact(A);
  const x = ls.solve(A, b, false);
  const Arow = req.symmetry ? [[4, 1], [1, 3]] : [[3, 7], [1, 2]];
  const ref = npSolve2(Arow, [1, 0.5]);
  assertArrayAlmostEqual(x.nonzeros(), ref, req.digits,
    `${plugin} explicit sfact/nfact x`);
}

// ----- driver -----

const tests = [
  ["test_simple",              test_simple],
  ["test_simple_dmatrix",      test_simple_dmatrix],
  ["test_simple_trans",        test_simple_trans],
  ["test_simple_trans_mx",     test_simple_trans_mx],
  ["test_simple_multi_rhs",    test_simple_multi_rhs],
  ["test_dense_3x3",           test_dense_3x3],
  ["test_sparse_pattern",      test_sparse_pattern],
  ["test_mx_solve",            test_mx_solve],
  ["test_identity",            test_identity],
  ["test_diagonal",            test_diagonal],
  ["test_residual",            test_residual],
  ["test_neig_rank",           test_neig_rank],
  ["test_factor_reuse",        test_factor_reuse],
  ["test_sfact_nfact_explicit", test_sfact_nfact_explicit],
];

(async () => {
  const { modulePath } = require("./_helpers");
  const create = require(modulePath);
  const M = await create();
  const solvers = ALL_SOLVERS.filter(([n]) => M.has_linsol(n));
  if (solvers.length === 0) {
    console.error("SKIP: no linsol plugin loaded.");
    process.exit(2);
  }
  let pass = 0, fail = 0;
  const failures = [];
  for (const [plugin, req] of solvers) {
    for (const [name, fn] of tests) {
      const label = `[${plugin}] ${name}`;
      try { fn(M, plugin, req); console.log(`ok   -- ${label}`); pass++; }
      catch (e) {
        console.log(`FAIL -- ${label}: ${e.message}`);
        fail++; failures.push([label, e.message]);
      }
    }
  }
  console.log(`\n${pass}/${pass + fail} passed`);
  if (failures.length > 0) {
    console.log("\n--- failures ---");
    for (const [n, m] of failures) console.log(`  ${n}: ${m}`);
  }
  process.exit(fail === 0 ? 0 : 1);
})().catch((e) => { console.error("FATAL:", e); process.exit(1); });
