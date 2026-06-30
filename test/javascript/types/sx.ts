// sx.ts -- Type-level coverage of the SX API surface.
//
// SX is the scalar-symbolic counterpart of MX (purely-symbolic
// expressions stored as a DAG of scalar SXElem operations).  Most of
// its surface mirrors MX, but a few methods (dep / n_dep / name /
// is_op / .elem(...)) reach into the per-scalar DAG and are SX-only.
// This fixture spot-checks the differences AND verifies the
// shared-surface paths (which would catch any SX-vs-MX divergence
// in the d.ts emitter).

import {
  SX, MX, DM, Sparsity,
  sin, cos, plus, times, minus, vertcat, mtimes,
  jacobian, gradient, substitute, dot, sumsqr,
} from "./casadi";

function expectType<T>(_v: T): void {}

// ============================================================
// Construction
// ============================================================
const s0 = SX();
const s1 = SX(3);
const sp = Sparsity.dense(2, 3);
const s2 = SX(sp);

expectType<SX>(SX.sym("x"));
expectType<SX>(SX.sym("x", 3));
expectType<SX>(SX.sym("x", 3, 2));
expectType<SX>(SX.sym("x", sp));
expectType<SX[]>(SX.sym("x", sp, 4));
expectType<SX[]>(SX.sym("x", 3, 2, 4));

// ============================================================
// Arithmetic / math
// ============================================================
const sx = SX.sym("sx", 3);
const sy = SX.sym("sy", 3);
expectType<SX>(plus(sx, sy));
expectType<SX>(plus(sx, 2));
expectType<SX>(times(sx, sy));
expectType<SX>(minus(sx, sy));
expectType<SX>(sin(sx));
expectType<SX>(cos(sx));
expectType<SX>(vertcat([sx, sy]));
expectType<SX>(mtimes(sx, sy.T()));

// ============================================================
// SX-only DAG introspection
// ============================================================
const expr = sin(sx);
expectType<bigint>(expr.n_dep());
expectType<SX>(expr.dep(0));
expectType<string>(expr.name());
expectType<boolean>(expr.is_op(0));

// ============================================================
// Differentiation
// ============================================================
const scalar = SX.sym("s");
const f_expr = sin(scalar);
expectType<SX>(jacobian(f_expr, scalar));
expectType<SX>(gradient(f_expr, scalar));
expectType<SX>(substitute(plus(sx, sy), sx, sy));
expectType<SX>(dot(sx, sy));
expectType<SX>(sumsqr(sx));

// ============================================================
// Cross-type: SX inputs accept DM-coercibles via _SX = SX | _DM
// ============================================================
const dx = DM([1, 2, 3]);
expectType<SX>(plus(sx, dx));  // _DM in _SX
expectType<SX>(plus(sx, 1.5));

// Pin locals
void [s0, s1, s2, sx, sy, expr, scalar, f_expr, dx];
