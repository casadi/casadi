// mx.ts -- Type-level coverage of the MX API surface used in
// test/python/mx.py.  Goal: any TS error that comes out of `tsc` here
// represents a real d.ts gap or clunkiness, NOT a fixture mistake.
// We do NOT add `// @ts-ignore` to silence things; we either fix the
// upstream typemap or accept the TS error as a guard against future
// regressions.

import {
  MX, DM, SX, Sparsity, Function as CFn,
  vertcat, horzcat, diagcat,
  mtimes, plus, minus, times, rdivide,
  sin, cos, tan, exp, log, sqrt, fmax, fmin, atan2,
  floor, ceil,
  reshape, repmat, transpose,
  jacobian, gradient, jtimes, hessian,
  substitute, dot,
} from "./casadi";

function expectType<T>(_v: T): void {}

// ============================================================
// MX construction
// ============================================================
const m0 = new MX();                       // empty
const m1 = new MX(0);                      // scalar 0
const sp = Sparsity.dense(3n, 2n);
const m2 = new MX(sp);                     // sparse pattern

// MX.sym in its full arity matrix
expectType<MX>(MX.sym("x"));
expectType<MX>(MX.sym("x", 3n));
expectType<MX>(MX.sym("x", 3n, 2n));
expectType<MX>(MX.sym("x", sp));
// Vector / matrix-of-MX symbolic primitives
expectType<MX[]>(MX.sym("x", sp, 4n));
expectType<MX[]>(MX.sym("x", 3n, 2n, 4n));
expectType<MX[][]>(MX.sym("x", sp, 4n, 5n));
expectType<MX[][]>(MX.sym("x", 3n, 2n, 4n, 5n));

// ============================================================
// Arithmetic (free-fn forms) -- coercion + return narrowing
// ============================================================
const x = MX.sym("x", 3n);
const y = MX.sym("y", 3n);

expectType<MX>(plus(x, y));
expectType<MX>(plus(x, 2));
expectType<MX>(plus(x, 2.5));
expectType<MX>(plus(2, x));
expectType<MX>(plus(x, [1, 2, 3]));
expectType<MX>(minus(x, y));
expectType<MX>(times(x, y));
expectType<MX>(rdivide(x, y));

// ============================================================
// Math functions
// ============================================================
expectType<MX>(sin(x));
expectType<MX>(cos(x));
expectType<MX>(tan(x));
expectType<MX>(exp(x));
expectType<MX>(log(x));
expectType<MX>(sqrt(x));
expectType<MX>(fmax(x, y));
expectType<MX>(fmin(x, y));
expectType<MX>(atan2(x, y));
expectType<MX>(floor(x));
expectType<MX>(ceil(x));

// ============================================================
// Concat / reshape
// ============================================================
expectType<MX>(vertcat([x, y]));
expectType<MX>(horzcat([x, y]));
expectType<MX>(diagcat([x, y]));
expectType<MX>(reshape(x, 1n, 3n));
expectType<MX>(reshape(x, sp));
expectType<MX>(repmat(x, 2n, 3n));
expectType<MX>(transpose(x));

// Method-form transpose (.T())
expectType<MX>(x.T());

// ============================================================
// Sparsity / shape introspection
// ============================================================
expectType<Sparsity>(x.sparsity());
expectType<bigint>(x.size1());
expectType<bigint>(x.size2());
expectType<bigint>(x.numel());
expectType<bigint>(x.nnz());
expectType<boolean>(x.is_dense());
expectType<boolean>(x.is_scalar());
expectType<boolean>(x.is_vector());
expectType<boolean>(x.is_square());

// ============================================================
// Function construction + evaluation
// ============================================================
const f = new CFn("f", [x, y], [plus(x, y), times(x, y)]);
expectType<CFn>(f);
expectType<bigint>(f.n_in());
expectType<bigint>(f.n_out());
expectType<Sparsity>(f.sparsity_in(0n));
expectType<Sparsity>(f.sparsity_in("i0"));

// Call: positional array form returns MX[]
expectType<MX[]>(f.call([x, y]));
// Call: dict form returns Record<string, MX>
expectType<Record<string, MX>>(f.call({ i0: x, i1: y }));

// Apply same Function to DM concrete values
const dx = DM.sym("dx", 3n);
const dy = DM.sym("dy", 3n);
expectType<DM[]>(f.call([dx, dy]));
expectType<Record<string, DM>>(f.call({ i0: dx, i1: dy }));

// ============================================================
// Differentiation
// ============================================================
const scalar = MX.sym("s");
const expr = sin(scalar);
expectType<MX>(jacobian(expr, scalar));
expectType<MX>(gradient(expr, scalar));
// hessian returns [H, g] tuple in C++; should reach via argout aggregation
const hg = hessian(expr, scalar);
expectType<[MX, MX]>(hg);
// jtimes: directional derivative
const v = MX.sym("v");
expectType<MX>(jtimes(expr, scalar, v));

// ============================================================
// Substitute / dot
// ============================================================
const a = MX.sym("a");
const b = MX.sym("b");
expectType<MX>(substitute(plus(a, b), a, b));
expectType<MX>(dot(x, y));

// ============================================================
// Indexing / slicing -- get / set with Slice + IM + sparsity
// ============================================================
import { Slice } from "./casadi";
const M = MX.sym("M", 3n, 3n);
// scalar element via two-index get
M.get(false, new Slice(0n, false), new Slice(1n, false));
// row slice
const Mt: MX = M.T();
expectType<MX>(Mt);

// ============================================================
// Serialization round-trip
// ============================================================
// MX has only `serialize(SerializingStream)` -- no parameterless form
// (C++ asymmetry: Sparsity / Matrix / Function DO have serialize() -> string
// but MX doesn't, see casadi/core/mx.hpp).  Use Sparsity instead.
const ser: string = sp.serialize();
expectType<string>(ser);

// ============================================================
// Function composition via call
// ============================================================
const f_unary = new CFn("f", [scalar], [sin(scalar)]);
const gx = f_unary.call([scalar]);
expectType<MX[]>(gx);

// Pin all locals so --noUnusedLocals stays quiet if ever enabled.
void [m0, m1, m2, M, Mt, ser, f_unary, gx, sp];
