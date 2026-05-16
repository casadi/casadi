// function.ts -- Type-level coverage of the Function API surface used
// in test/python/function.py.  Goal: surface d.ts gaps via real
// consumer patterns; no @ts-ignore.

import {
  Function as CFn, MX, SX, DM, Sparsity,
  vertcat, plus, times, sin,
} from "./casadi";

function expectType<T>(_v: T): void {}

// ============================================================
// Construction
// ============================================================
const x = MX.sym("x");
const y = MX.sym("y");
const f = new CFn("f", [x, y], [plus(x, y), times(x, y)]);
expectType<CFn>(f);

// With named input/output
const f2 = new CFn("f2", [x, y], [plus(x, y)], ["a", "b"], ["s"]);
expectType<CFn>(f2);

// SX version
const sx0 = SX.sym("sx");
const sf = new CFn("sf", [sx0], [sin(sx0)]);
expectType<CFn>(sf);

// JIT-compiled scalar function from source string
const jf = CFn.jit("jf", "function jf(x: number): number { return 2*x; }", ["x"], ["y"]);
expectType<CFn>(jf);

// ============================================================
// Introspection
// ============================================================
expectType<bigint>(f.n_in());
expectType<bigint>(f.n_out());
expectType<string[]>(f.name_in());
expectType<string[]>(f.name_out());
expectType<string>(f.name_in(0n));
expectType<string>(f.name_out(0n));
expectType<bigint>(f.index_in("a"));
expectType<bigint>(f.index_out("s"));
expectType<Sparsity>(f.sparsity_in(0n));
expectType<Sparsity>(f.sparsity_in("a"));
expectType<Sparsity>(f.sparsity_out(0n));
expectType<[bigint, bigint]>(f.size_in(0n));
expectType<[bigint, bigint]>(f.size_in("a"));
expectType<[bigint, bigint]>(f.size_out(0n));
expectType<bigint>(f.nnz_in());
expectType<bigint>(f.nnz_in(0n));
expectType<bigint>(f.numel_in());
expectType<bigint>(f.numel_in(0n));
// info() returns a primitive-valued dict (GenericType variants are
// unwrapped at the JS marshaling boundary; see %ts_alias_out wiring).
// Verify the returned-type allows JSON-shaped destructuring.
const info = f.info();
const v = info["some_key"];
// v is the union of primitive types; verify we can narrow:
if (typeof v === "number") expectType<number>(v);
if (typeof v === "string") expectType<string>(v);
if (typeof v === "boolean") expectType<boolean>(v);
if (typeof v === "bigint") expectType<bigint>(v);

// Symbolic input/output primitive enumeration
expectType<MX>(f.mx_in(0n));
expectType<MX[]>(f.mx_in());
expectType<SX>(sf.sx_in(0n));
expectType<SX[]>(sf.sx_in());

// ============================================================
// Calling
// ============================================================
// Positional list-of-MX
const out_pos = f.call([x, y]);
expectType<MX[]>(out_pos);
// Dict form with MX
const out_dict = f.call({ a: x, b: y });
expectType<Record<string, MX>>(out_dict);
// Numerical evaluation with DM
const dx = new DM([1.0]);
const dy = new DM([2.0]);
const out_num = f.call([dx, dy]);
expectType<DM[]>(out_num);
// Numerical evaluation with dict
const out_dict_num = f.call({ a: dx, b: dy });
expectType<Record<string, DM>>(out_dict_num);

// always_inline / never_inline flags
expectType<MX[]>(f.call([x, y], true));
expectType<MX[]>(f.call([x, y], false, true));

// ============================================================
// Derivative builders
// ============================================================
expectType<CFn>(f.jacobian());
expectType<CFn>(f.expand());
expectType<CFn>(f.expand("f_sx"));

// ============================================================
// mapaccum / fold / map
// ============================================================
expectType<CFn>(f.mapaccum(10n));
expectType<CFn>(f.mapaccum("f_acc", 10n));
expectType<CFn>(f.fold(10n));
expectType<CFn>(f.map(10n));
expectType<CFn>(f.map(10n, "openmp"));
expectType<CFn>(f.map(10n, "thread", 4n));

// ============================================================
// Codegen
// ============================================================
expectType<string>(f.generate("f_codegen"));
expectType<string>(f.generate("f_codegen", { with_header: true }));

// ============================================================
// Stats (post-eval)
// ============================================================
// stats() has the same shape as info() -- primitive-valued dict.
const s = f.stats();
const iter_count = s["iter_count"];
if (typeof iter_count === "bigint") expectType<bigint>(iter_count);

// Pin locals
void [f, f2, sf, jf, out_pos, out_dict, out_num, out_dict_num, info, v, s, iter_count];
