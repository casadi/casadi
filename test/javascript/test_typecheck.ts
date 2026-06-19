// test_typecheck.ts
//
// TYPE-ONLY smoke test for the auto-generated casadi.d.ts emitted by
// `swig -wasm-js -stubs`.  We only exercise the type system: tsc
// --noEmit must accept this file without error, but nothing runs at
// runtime.  The wasm module is loaded by `test_runtime.cjs`, not here.
//
// Every assignment below pins a concrete TypeScript type against a
// method's declared return type -- a regression in the tsstub_in /
// tsstub_out plumbing (or in cpp_to_ts fallbacks) flips the matching
// assignment into a tsc error.

import type {
  Function as Func,
  SX,
  MX,
  DM,
  Sparsity,
  Slice,
  GenericType,
  Dict,
} from "../../build-wasm/swig/wasm-js/casadi";

// ---- Class instances (declared, not constructed) ------------------------
declare const f: Func;
declare const sx: SX;
declare const mx: MX;
declare const dm: DM;

// ---- Primitive return-type pinning --------------------------------------
const n_in:   bigint  = f.n_in();
const n_out:  bigint  = f.n_out();
const fname:  string  = f.name();

// ---- Container return types ---------------------------------------------
const ins:  string[] = f.name_in();
const outs: string[] = f.name_out();

// ---- Methods returning class instances ----------------------------------
const expanded: Func = f.expand();

// ---- Dict type alias (from casadi.i %insert("stubs") preamble) ----------
declare const d: Dict;
const dict: Record<string, GenericType> = d;

// ---- Negative checks ----------------------------------------------------
// Each `@ts-expect-error` MUST produce a tsc error; if any compiles
// the assertion below it would silently pass invalid types.

// @ts-expect-error -- n_in returns bigint, not number
const bad_nin: number = f.n_in();

// @ts-expect-error -- name() returns string, not bigint
const bad_name: bigint = f.name();

// @ts-expect-error -- expand() returns Function, not string
const bad_expand: string = f.expand();

// @ts-expect-error -- name_in() returns string[], not bigint[]
const bad_ins: bigint[] = f.name_in();

void [n_in, n_out, fname, ins, outs, expanded, dict,
      bad_nin, bad_name, bad_expand, bad_ins, sx, mx, dm];
