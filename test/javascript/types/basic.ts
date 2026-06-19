// basic.ts -- Sanity fixture for the SWIG-generated casadi.d.ts.
//
// Exercises the core happy-path patterns that callers should be able
// to write without any manual wrapping or type-assertions:
//
//   * Typed static factories on the leaf class (e.g. MX.sym not
//     MatrixCommon.sym).
//   * Return-type narrowing (sin(DM) -> DM, not any).
//   * Constructors with their full overload set.
//   * Interface-merged base methods reachable on the subclass.
//
// `expectType<T>(value)` is the TypeScript analog of the
// `typing_extensions.assert_type` calls in test/python/pyright_stubs.py:
// any T-vs-inferred-type mismatch surfaces as a tsc error at the call
// site, which the tsc_stubs.js harness counts against the budget.

import {
  DM, SX, MX, Sparsity, Linsol,
  MatrixCommon, GenericExpressionCommon, PrintableCommon,
} from "./casadi";

function expectType<T>(_value: T): void { /* type-only assertion */ }

// --- Constructors ---
const d0 = DM();
const d1 = DM(5);
const d2 = DM([[1, 2], [3, 4]]);
const sp = Sparsity.dense(3, 3);
const d3 = DM(sp);

// --- Static `sym` reachable on the subclass (was previously only on Gen<X>) ---
const dx = DM.sym("dx", 3);
const sx = SX.sym("sx", 3);
const mx = MX.sym("mx", 3);
expectType<DM>(dx);
expectType<SX>(sx);
expectType<MX>(mx);

// --- Sparsity constructor variant ---
const sx2 = SX(sp);
expectType<SX>(sx2);

// --- Linsol ---
const ls = Linsol("S", "qr", sp);
expectType<Linsol>(ls);

// --- Interface-merged base reachable through declaration merging ---
const _base_dm: MatrixCommon = dx;
const _base_mx: GenericExpressionCommon = mx;
const _print_mx: PrintableCommon = mx;

// Keep references so unused-locals-strict doesn't complain.
void [d0, d1, d2, d3, sx2, ls, _base_dm, _base_mx, _print_mx];
