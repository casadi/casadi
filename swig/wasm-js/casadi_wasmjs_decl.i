// Linsol: only as the element of a std::vector in the serializer API
%fragment(SWIG_Traits_frag(casadi::Linsol));

// Director typemaps; after %casadi_typemaps (last registered wins)

%casadi_director_class(casadi::Sparsity)
%casadi_director_class(casadi::Matrix<double>)
%casadi_director_class(casadi::MX)
%casadi_director_class(casadi::Matrix<casadi::SXElem>)
%casadi_director_class(casadi::Matrix<casadi_int>)
%casadi_director_class(casadi::Function)
%casadi_director_class(casadi::Slice)
%casadi_director_class(casadi::SXElem)
// No GenericType: JS sees plain values

%casadi_director_vec(casadi::Sparsity)
%casadi_director_vec(casadi::Matrix<double>)
%casadi_director_vec(casadi::MX)
%casadi_director_vec(casadi::Matrix<casadi::SXElem>)
%casadi_director_vec(casadi::Function)

%casadi_director_vec_str

// jsout for every dict-shaped return type
%casadi_dict_jsout(casadi::GenericType)
%casadi_dict_jsout(casadi::Sparsity)
%casadi_dict_jsout(casadi::Matrix<double>)
%casadi_dict_jsout(casadi::MX)
%casadi_dict_jsout(casadi::Matrix<casadi::SXElem>)
%casadi_dict_jsout(std::vector<std::string>)
%casadi_dict_jsout(std::vector< casadi::Sparsity >)
%casadi_dict_jsout(std::vector< casadi::Matrix<double> >)
%casadi_dict_jsout(std::vector< casadi::MX >)
%casadi_dict_jsout(std::vector< casadi::Matrix<casadi::SXElem> >)

%insert("stubs") %{
// --- casadi-specific TypeScript type aliases ---
export type Dict         = Record<string, GenericType>;
export type DMDict       = Record<string, DM>;
export type MXDict       = Record<string, MX>;
export type SXDict       = Record<string, SX>;
export type SparsityDict = Record<string, Sparsity>;

// --- Input-position coercion unions (analog of Python pyi's _DM / _SX / _MX).
// xTsStub args in casadi.i use the bare class names (DM / SX / MX / ...);
// the wasm_js.cxx d.ts emitter rewrites bare identifiers to their `_`
// aliases at input positions only, so users can pass JS-natural values
// (numbers, arrays, Sparsity patterns) where DM/SX/MX is expected.  Output
// positions stay as the bare class (DM/SX/MX) so return-type narrowing
// works.
export type _DM = boolean | bigint | number | DM | Sparsity | number[] | number[][];
export type _SX = SX | _DM;
export type _MX = MX | _DM;
export type _IM = IM | bigint | number | number[] | number[][];
export type _SXElem = SXElem | number;
export type _Slice = Slice | number | bigint;
// TS 4.1+ permits inline array/object types in recursive aliases, but a
// bare reference like `_GenericType[]` triggers TS2456 (circular).  Use
// the inline-array form to break the cycle the compiler can see through.
export type _GenericType = boolean | bigint | number | string | Function
                         | readonly _GenericType[]
                         | { readonly [k: string]: _GenericType };
// JS callers regularly pass plain `number` literals where casadi_int
// (-> bigint on the wasm boundary) is expected.  The wasm marshaling
// auto-converts, so widen `bigint` -> `_bigint = bigint | number` on
// INPUT positions only.  Output stays `bigint` for narrow typing of
// size / index / count returns.  Hooked up via %ts_alias_in(bigint)
// below.
export type _bigint = bigint | number;
export type _MIndex = number | bigint | bigint[] | boolean[] | Sparsity | DM;

// Note: DMVector / MXVector / SXVector are NOT type aliases here -- they
// collide with the SWIG-emitted std::vector wrapper classes at the top
// of the d.ts (export class DMVector { ... }).  Functions accepting
// std::vector<DM> use `_DM[]` (input) / `DM[]` (output) via the emitter's
// coercion-union rewrite.  Users who want the wrapper class still get
// `new DMVector(arr)` via the class declaration.
export type DMVectorDict = Record<string, DM[]>;
export type MXVectorDict = Record<string, MX[]>;
export type SXVectorDict = Record<string, SX[]>;
export type IM           = any;  // casadi::Matrix<casadi_int>: no JS-side class today
export type casadi_index = [bigint, bigint];
export type native_DM    = Float64Array;
// Serialization + internal-impl classes -- not exposed to JS.
export type SerializingStream    = any;
export type DeserializingStream  = any;
export type SharedObjectInternal = any;
// Enums that aren't auto-wrapped.
export type ConstraintType = number;
export type VariableType   = number;
export type DomainType     = number;

// ---------------------------------------------------------------------
// JS-side concat aliases.  At runtime, casadi.i's %insert("js")
// block adds variadic + list-form aliases to the concat family:
//   vertcat(a, b, c)  -- variadic
//   vcat([a, b, c])   -- list form
//   (analogous hcat / dcat).
// The auto-generated d.ts only declares the C++-direct names
// (vertcat / horzcat / diagcat).  Re-export the list-form aliases
// so JS-test code calling `M.vcat([...])` typechecks.
export declare const vcat: typeof vertcat;
export declare const hcat: typeof horzcat;
export declare const dcat: typeof diagcat;

// Make Function instances callable at the type level.  At runtime,
// casadi.js wraps Function instances so `f(arg)` evaluates as
// shorthand for `f.call(arg)`.  The auto-generated class body only
// has `.call(...)` so without this merge, `f([3, 7])` would
// typecheck as "instance has no call signature".  This is a
// casadi-specific runtime convention that the wasm_js.cxx emitter
// has no business knowing about.
export interface Function__class {
  (arg: _DM[]): DM[];
  (arg: _SX[]): SX[];
  (arg: _MX[]): MX[];
  (arg: Record<string, _DM>): Record<string, DM>;
  (arg: Record<string, _SX>): Record<string, SX>;
  (arg: Record<string, _MX>): Record<string, MX>;
  // Positional varargs (the runtime accepts `f(x, y, z)` as
  // shorthand for `f([x, y, z])`).  Loose return type since the
  // input flavor is heterogeneous in this form.
  (...args: any[]): any;
}

// Runtime factory shape.
//
// At runtime, `require("./casadi")` returns an async factory function
// (the SWIG `module.exports = async function createcasadi() {...}`
// pattern, see build-wasm/swig/wasm-js/casadi.js).  The factory loads
// the underlying wasm module and returns an object whose properties
// are ALL the exports declared above (`DM`, `SX`, `vertcat`, ...).
//
// Static .ts consumers use the named-export form
// (`import { DM } from "./casadi"`) and never call the factory --
// that path is fine because the d.ts shape AND the JS runtime shape
// happen to agree on the *names* of the exports.
//
// Runtime .js consumers (and TS callers that want to invoke the
// factory) use:
//   const create = require("./casadi");
//   const ca = await create();
//   ca.DM.sym("x");
// To make `ca` typed against the same export shape, the factory's
// return-type is declared self-referentially via
// `typeof import("./casadi")`.  No double bookkeeping needed.
declare const createcasadi: () => Promise<typeof import("./casadi")>;
export default createcasadi;
%}

// d.ts: inputs widen to _<name>, outputs stay bare (GenericType unwraps)
%ts_alias_in(DM)
%ts_alias_in(SX)
%ts_alias_in(MX)
%ts_alias_in(IM)
%ts_alias_in(SXElem)
%ts_alias_in(Slice)
%ts_alias_in(GenericType)
%ts_alias_in(MIndex)
%ts_alias_in(bigint)

%ts_alias_out(GenericType)

// JS-side ergonomics, emitted after the class definitions
%insert("js") %{
  /* Coerce bare numbers/arrays into matrix proxies (jsarg typemaps),
     so `plus(x, 2)` works. */
  const __coerce = (v, Cls) =>
    (v === null || v === undefined ||
     (typeof v === "object" && v._swig_type !== undefined)) ? v : new Cls(v);

  /* Element indexing (wired by %feature("wasmjs:index")), mirroring
     Python: x[k], x[k]=v, x["a:b"], x["0::3"], x["a:b,j"], x.nz[k]
     all route to casadi's 0-based get/set/get_nz/set_nz. */
  const __isIdxKey = (k) => typeof k === "string" && k.length > 0 &&
    (k.charCodeAt(0) === 45 || k.charCodeAt(0) === 58 ||  // '-' or ':'
     (k.charCodeAt(0) >= 48 && k.charCodeAt(0) <= 57));   // digit
  const __idxPart = (s) => {
    if (s.indexOf(":") < 0) return Number(s);
    if (s === ":") return new Slice();
    const p = s.split(":");
    /* Python-style open bounds: "0::3" / ":b" / "1:".  Empty start -> 0;
       empty stop -> the Slice unboundedness sentinel (int64 max, see
       slice.hpp) -- exactly what Python maps a `None` stop to. */
    const start = p[0] === "" ? 0 : Number(p[0]);
    const stop = p[1] === "" ? 9223372036854775807n : Number(p[1]);
    return p.length >= 3 ? new Slice(start, stop, Number(p[2]))
                         : new Slice(start, stop);
  };
  const __idxArgs = (k) => k.split(",").map(__idxPart);
  const __NZ_HANDLER = {
    get(t, k, r) { return __isIdxKey(k) ? t.get_nz(false, Number(k)) : Reflect.get(t, k, r); },
    set(t, k, v, r) { if (__isIdxKey(k)) { t.set_nz(v, false, Number(k)); return true; } return Reflect.set(t, k, v, r); },
  };
  const __INDEX_HANDLER = {
    get(t, k, r) {
      if (k === "__indexed") return true;
      if (k === "nz") return new Proxy(t, __NZ_HANDLER);
      if (__isIdxKey(k)) return t.get(false, ...__idxArgs(k));
      return Reflect.get(t, k, r);
    },
    set(t, k, v, r) {
      if (__isIdxKey(k)) { t.set(v, false, ...__idxArgs(k)); return true; }
      return Reflect.set(t, k, v, r);
    },
  };
  const __mkIndex = (o) =>
    (o === null || typeof o !== "object" || o.__indexed) ? o : new Proxy(o, __INDEX_HANDLER);

  /* x.T as a getter, matching Python. */
  for (const Cls of [DM, MX, SX, Sparsity]) {
    const __orig_T = Cls.prototype.T;
    if (!__orig_T) continue;
    Object.defineProperty(Cls.prototype, "T", {
      get() { return __orig_T.call(this); },
      configurable: true,
    });
  }

  /* Concat family, mirroring Python: vertcat(a, b, c) is variadic,
     vcat([a, b, c]) takes a list (same for horzcat/hcat, diagcat/dcat). */
  const __cat_variadic = (orig) => function (...args) {
    return orig.call(this, args);
  };
  const __cat_list = (orig) => function (arr) {
    return orig.call(this, arr);
  };
  for (const [variadic, list] of [
    ["horzcat", "hcat"],
    ["vertcat", "vcat"],
    ["diagcat", "dcat"],
  ]) {
    const fn = __m[variadic];
    __m[variadic] = __cat_variadic(fn).bind(__m);
    __m[list]     = __cat_list(fn).bind(__m);
  }

  /* Mixed input lists ([DM(0), mx, ...]) must dispatch as the widest
     element type; the generated per-element probes are exact-class, so
     pre-coerce JS-side (matches Python's implicit DM->MX/SX promotion). */
  const __origFnCall = Function.prototype.call;
  Function.prototype.call = function (arg, ...rest) {
    if (Array.isArray(arg)) {
      if (arg.some((e) => e instanceof MX)) arg = arg.map((e) => e instanceof MX ? e : new MX(e));
      else if (arg.some((e) => e instanceof SX)) arg = arg.map((e) => e instanceof SX ? e : new SX(e));
    }
    return __origFnCall.call(this, arg, ...rest);
  };

  // Preserve symbolic lists before permissive scalar-matrix overloads accept arrays.
  const __origExportGraph = __m.export_graph;
  __m.export_graph = function (expression, ...args) {
    if (Array.isArray(expression)) {
      const Cls = expression.every(e => e instanceof SX) ? SXVector
        : expression.every(e => e instanceof MX) ? MXVector : null;
      if (Cls) {
        const vector = __arr_to_vec(expression, Cls);
        try { return __origExportGraph.call(this, vector, ...args); }
        finally { vector.delete(); }
      }
    }
    return __origExportGraph.call(this, expression, ...args);
  };

  /* Function instances are directly callable, mirroring Python:
     f(a, b) positional with n_out-shaped return (0 -> null, 1 -> bare,
     >1 -> list); f({name: v}) dict-in/dict-out.  The dict check must
     exclude casadi types, which are objects too.  Generic machinery
     (prototype transfer, FinalizationRegistry handoff, re-wrap of
     methods returning bare instances) is __callable_class, backend-emitted. */
  Function = __m.Function = __callable_class(Function, __fr_casadi_Function, (f) => function (...args) {
    if (args.length === 1) {
      const a = args[0];
      if (a && typeof a === "object" && !Array.isArray(a)
          && !(a instanceof DM) && !(a instanceof MX) && !(a instanceof SX)
          && !(a instanceof DMVector) && !(a instanceof MXVector) && !(a instanceof SXVector)
          && !(a instanceof Sparsity) && (typeof IM === "undefined" || !(a instanceof IM)))
        return f.call(a);
    }
    const r = f.call(args);
    const n = Number(f.n_out());
    return n === 0 ? null : n === 1 ? r[0] : r;
  });

%}
%insert("js") "wasm-js/load_plugins.js"


// Enum values as 64-bit ints (bigint in JS, like casadi_int)
%feature("wasmjs:enum_int_type","long long");

// Accept a bare JS number/array where SX/MX/DM is expected
%typemap(jstypecheck) casadi::Matrix<casadi::SXElem>, const casadi::Matrix<casadi::SXElem>&, casadi::Matrix<casadi::SXElem>&
  "(M._swig_can_SX(__unwrap($input)) || typeof $input === 'number' || typeof $input === 'bigint' || Array.isArray($input))";
%typemap(jsarg) casadi::Matrix<casadi::SXElem>, const casadi::Matrix<casadi::SXElem>&, casadi::Matrix<casadi::SXElem>&
  "__unwrap(__coerce($input, SX))";
%typemap(jstypecheck) casadi::MX, const casadi::MX&, casadi::MX&
  "(M._swig_can_MX(__unwrap($input)) || typeof $input === 'number' || typeof $input === 'bigint' || Array.isArray($input))";
%typemap(jsarg) casadi::MX, const casadi::MX&, casadi::MX&
  "__unwrap(__coerce($input, MX))";
%typemap(jstypecheck) casadi::Matrix<double>, const casadi::Matrix<double>&, casadi::Matrix<double>&
  "(M._swig_can_DM(__unwrap($input)) || typeof $input === 'number' || typeof $input === 'bigint' || Array.isArray($input))";
%typemap(jsarg) casadi::Matrix<double>, const casadi::Matrix<double>&, casadi::Matrix<double>&
  "__unwrap(__coerce($input, DM))";

// Element indexing via the __mkIndex Proxy
%feature("wasmjs:index") casadi::Matrix<casadi::SXElem>;
%feature("wasmjs:index") casadi::MX;
%feature("wasmjs:index") casadi::Matrix<double>;

// Opti values return as DM
%typemap(out, doc="DM", tsstub_out="DM", noblock=1, fragment="casadi_all") casadi::native_DM {
  if(!($result = swig::from< casadi::DM >($1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type 'DM'.");
}
