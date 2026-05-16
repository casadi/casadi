// typemaps.ts -- Type-level coverage of the SWIG typemap surface from
// the JS side.  Inspired by test/python/typemaps.py, which exercises
// numpy/scipy interop + Python operator-overloading autoconversion.
// JS doesn't have those, but the underlying typemap library still
// drives:
//   * DM/SX/MX construction from primitives (number, number[], number[][])
//   * Cross-type construction (SX from DM, MX from DM)
//   * Sparsity construction from row/col triplets, dense, banded
//   * Function evaluation with mixed input types
//   * INOUT typemaps (where they exist)
// Goal: surface any remaining d.ts edge case in the input-coercion
// unions (_DM / _SX / _MX) by exercising every C++ typemap branch.

import {
  DM, SX, MX, Sparsity, Function as CFn,
  plus, minus, times, rdivide, mtimes,
  vertcat, horzcat, blockcat,
  jacobian, gradient,
} from "./casadi";

function expectType<T>(_v: T): void {}

// ============================================================
// DM construction -- every typemap branch the C++ exposes
// ============================================================
expectType<DM>(new DM());                                // empty
expectType<DM>(new DM(5));                               // scalar
expectType<DM>(new DM(2.3));                             // scalar (float)
expectType<DM>(new DM([1, 2, 3]));                       // 1D number[]
expectType<DM>(new DM([[1, 2], [3, 4]]));                // 2D number[][]
expectType<DM>(new DM(Sparsity.dense(3n, 2n)));          // sparsity-only
expectType<DM>(new DM(Sparsity.dense(3n, 2n), new DM(1)));  // (sp, val)
expectType<DM>(new DM(new DM(7)));                       // copy

// ============================================================
// DM nonzeros + readback
// ============================================================
const d = new DM([[1, 2, 3], [4, 5, 6]]);
expectType<number[]>(d.nonzeros());
expectType<[bigint, bigint]>(d.size());
expectType<bigint>(d.size1());
expectType<bigint>(d.size2());
expectType<bigint>(d.numel());
expectType<bigint>(d.nnz());

// ============================================================
// Cross-type construction (SX from DM, MX from DM, etc.)
// ============================================================
expectType<SX>(new SX(new DM([1, 2, 3])));
expectType<MX>(new MX(new DM([1, 2, 3])));
expectType<MX>(new MX(5));        // scalar
expectType<MX>(new MX(2.5));      // float
// Note: SX(number[]) and MX(number[]) require Sparsity arg in C++ --
// the constructor surface is asymmetric.  Document the workaround:
expectType<SX>(new SX(new DM([1, 2, 3])));   // round-trip via DM
expectType<MX>(new MX(new DM([1, 2, 3])));

// ============================================================
// Autoconversion via _DM coercion union (input position widening)
// ============================================================
const x = MX.sym("x");
expectType<MX>(plus(x, 2.3));                            // number coerced
expectType<MX>(plus(x, [1, 2, 3]));                      // number[] coerced
expectType<MX>(plus(x, [[1, 2], [3, 4]]));               // number[][] coerced
expectType<MX>(plus(x, Sparsity.dense(1n)));             // Sparsity coerced
expectType<MX>(times(2.3, x));                           // commute
expectType<MX>(times([1, 2], x));
expectType<MX>(minus(x, 1));
expectType<MX>(rdivide(x, 2));
expectType<MX>(rdivide(2, x));

// ============================================================
// mtimes -- matrix multiplication, both 2-arg and array forms
// ============================================================
const A = new DM([[1, 2], [3, 4]]);
const B = new DM([[5, 6], [7, 8]]);
expectType<DM>(mtimes(A, B));                            // 2-arg
expectType<DM>(mtimes([A, B]));                          // array form
expectType<DM>(mtimes([A, B, A]));                       // chain

// With number coercion
expectType<DM>(mtimes(A, 2));

// ============================================================
// Sparsity construction from triplets
// ============================================================
expectType<Sparsity>(Sparsity.triplet(3n, 3n, [0n, 1n], [0n, 2n]));
// triplet with invert_mapping returns tuple
const [sp_tri, mapping] = Sparsity.triplet(3n, 3n, [0n, 1n], [0n, 2n], true);
expectType<Sparsity>(sp_tri);
expectType<bigint[]>(mapping);

// DM/SX triplet (values-keyed)
expectType<DM>(DM.triplet([0n, 1n], [0n, 2n], new DM([1.0, 2.0])));
expectType<DM>(DM.triplet([0n, 1n], [0n, 2n], new DM([1.0, 2.0]), 3n, 3n));
expectType<SX>(SX.triplet([0n, 1n], [0n, 2n], SX.sym("v", 2n)));

// ============================================================
// blockcat -- both varargs-4 and array-2D forms
// ============================================================
expectType<DM>(blockcat(A, B, A, B));                    // 4-arg form
expectType<DM>(blockcat([[A, B], [A, B]]));              // 2D array form

// ============================================================
// vertcat / horzcat across types
// ============================================================
expectType<DM>(vertcat([A, A]));
expectType<DM>(vertcat([A, A]));
expectType<MX>(vertcat([x, x]));
expectType<Sparsity>(vertcat([Sparsity.dense(2n), Sparsity.dense(2n)]));

// ============================================================
// Function eval with mixed-type inputs
// ============================================================
const u = SX.sym("u");
const f = new CFn("f", [u], [plus(u, new SX(1))]);
// Symbolic call: SX[] in -> SX[] out
const r_sym = f.call([u]);
expectType<SX[]>(r_sym);
// Numerical call: pass DM, get DM
const r_num = f.call([new DM(2)]);
expectType<DM[]>(r_num);
// Dict-form: same flavor in/out
const r_dict = f.call({ i0: new DM(2) });
expectType<Record<string, DM>>(r_dict);

// ============================================================
// Function eval with chained MX expressions
// ============================================================
const xMX = MX.sym("xMX");
const fMX = new CFn("fMX", [xMX], [plus(xMX, new MX(1))]);
const r_chain = fMX.call([plus(xMX, xMX)]);  // pass MX expr as input
expectType<MX[]>(r_chain);

// ============================================================
// Differentiation typemaps -- expr + var return correct flavor
// ============================================================
expectType<MX>(jacobian(plus(xMX, xMX), xMX));
expectType<MX>(gradient(plus(xMX, xMX), xMX));
expectType<SX>(jacobian(plus(u, u), u));
expectType<SX>(gradient(plus(u, u), u));

// ============================================================
// DM-SX-MX cast via Function eval (test_DMSX / test_DMMX patterns)
// ============================================================
const w_dm = new DM([[1, 2, 3], [4, 5, 6]]);

// Symbolic input picks up the constant via SX flavor
const x_sx = SX.sym("x");
const f_const_sx = new CFn("f", [x_sx], [new SX(w_dm)]);
// f.sx_in() returns SX[] (the symbolic inputs)
const sxin: SX[] = f_const_sx.sx_in();
expectType<SX[]>(sxin);
// f.call(f.sx_in()) returns SX[] with the constant lifted to SX flavor
const w_sx_arr = f_const_sx.call(sxin);
expectType<SX[]>(w_sx_arr);
expectType<bigint>(w_sx_arr[0].size1());

// Same with MX flavor
const x_mx = MX.sym("x");
const f_const_mx = new CFn("f", [x_mx], [new MX(w_dm)]);
const mxin: MX[] = f_const_mx.mx_in();
expectType<MX[]>(mxin);
const w_mx_arr = f_const_mx.call(mxin);
expectType<MX[]>(w_mx_arr);

// ============================================================
// set / get / set_nz / get_nz (DM mutation surface)
// ============================================================
const target = new DM([[0, 0], [0, 0]]);
target.set(new DM(5), false, 0n, 0n);   // set scalar at (0,0)
target.set(new DM([1, 2]), false, 0n);  // set row
expectType<DM>(target.get(false, 0n, 0n));
expectType<DM>(target.get(false, 0n));
expectType<DM>(target.get_nz(false, 0n));
target.set_nz(new DM(7), false, 0n);

void [d, x, A, B, sp_tri, mapping, u, f, r_sym, r_num, r_dict, xMX, fMX, r_chain,
      w_dm, x_sx, f_const_sx, sxin, w_sx_arr, x_mx, f_const_mx, mxin, w_mx_arr, target];
