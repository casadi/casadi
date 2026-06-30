// overloads.ts -- Fixture for free-function overload resolution.
//
// The d.ts emitter buckets free-function overloads per name and sorts
// each bucket by typecheck precedence ascending so TS picks the
// most-specific matching overload.  Without the sort, `plus(dm, 2)`
// would match the `plus(_MX, _MX): MX` overload first (because
// `_MX = MX | _DM` includes DM-coercible types) and silently return
// `MX` instead of `DM`.  These expectType checks regress on the wrong
// return type.

import { DM, SX, MX, Sparsity, plus, mtimes, vertcat, horzcat, sin, cos } from "./casadi";

function expectType<T>(_value: T): void {}

const dx = DM.sym("dx", 3);
const sx = SX.sym("sx", 3);
const mx = MX.sym("mx", 3);
const sp = Sparsity.dense(3, 3);

// --- Most-specific overload selected per leading-arg type ---
expectType<DM>(plus(dx, dx));
expectType<SX>(plus(sx, sx));
expectType<MX>(plus(mx, mx));
// `plus` is element-wise math; no Sparsity-only overload exists.
// Use vertcat / horzcat for the Sparsity-overload exercise.
expectType<Sparsity>(vertcat([sp, sp]));
expectType<Sparsity>(horzcat([sp, sp]));

// --- Input coercion: number / number[] / Sparsity widened via _DM ---
expectType<DM>(plus(dx, 2));
expectType<DM>(plus(dx, 3.5));
expectType<DM>(plus(dx, [1, 2, 3]));
expectType<DM>(plus(dx, sp));
expectType<DM>(mtimes(dx, 4));

// --- _SX / _MX widening includes _DM ---
expectType<SX>(plus(sx, 1));
expectType<MX>(plus(mx, 1));

// --- Heterogeneous: DM + MX picks MX (DM is in _MX) ---
expectType<MX>(plus(mx, dx));
expectType<MX>(mtimes(mx, dx));

// --- Unary functions narrow correctly ---
expectType<DM>(sin(dx));
expectType<SX>(sin(sx));
expectType<MX>(sin(mx));
expectType<number>(sin(0.5));
expectType<number>(cos(1.0));

// --- Array overloads ---
expectType<DM>(vertcat([dx, dx]));
expectType<SX>(vertcat([sx, sx]));
expectType<MX>(vertcat([mx, mx]));
expectType<DM>(horzcat([dx, dx]));
