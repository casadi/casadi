// Smoke test for the newly-bridged Callback director methods (post-Layer 1-3
// refactor): get_jacobian, get_forward, get_reverse take const Dict& opts +
// const std::vector<std::string>& names which previously had no wasm-js
// directorin and got Warning 469 / skipped at SWIG-time.
const path = require("path");
const create = require(path.resolve(
  "/home/jgillis/programs/casadi/build-wasm-flang/build/swig/wasm-js/casadi.js"
));

(async () => {
  const M = await create();
  let pass = 0, fail = 0;

  // --- 1. get_jacobian override receives (name, inames, onames, opts) ---
  let jac_call_log = null;
  class JacCB extends M.Callback {
    constructor() {
      super();
      this.construct("JacCB", null);
    }
    get_n_in()  { return 1n; }
    get_n_out() { return 1n; }
    has_jacobian() { return true; }
    get_jacobian(args) {
      // Director calling convention (matches existing eval(args) in
      // director.js): args is a single JS array containing the marshalled
      // C++ parms in order -- here (name, inames, onames, opts).
      const name = args[0], inames = args[1], onames = args[2], opts = args[3];
      jac_call_log = {
        name_type:    typeof name,
        name:         name,
        inames_arr:   Array.isArray(inames),
        inames_len:   inames.length,
        onames_arr:   Array.isArray(onames),
        opts_type:    typeof opts,
      };
      // Return a Function that just computes 2*x (Jacobian of x^2 mock).
      const x = M.SX.sym("x_in");
      const out = M.SX.sym("nominal_out");
      const jac = M.times(M.SX(2), x);
      return M.Function(name, [x, out], [jac]);
    }
    eval(args) {
      const x = args[0].nonzeros()[0];
      return [M.DM(x * x)];
    }
  }
  const cb = new JacCB();
  try {
    // Callback is-a Function, so jacobian() is directly callable -- it
    // drives Function::factory -> get_jacobian under the hood when
    // has_jacobian() returns true.
    cb.jacobian();
    if (!jac_call_log) {
      console.log("FAIL get_jacobian: never invoked");
      fail++;
    } else if (jac_call_log.name_type !== "string") {
      console.log("FAIL get_jacobian: name not a string:", jac_call_log);
      fail++;
    } else if (!jac_call_log.inames_arr) {
      console.log("FAIL get_jacobian: inames not an Array:", jac_call_log);
      fail++;
    } else if (!jac_call_log.onames_arr) {
      console.log("FAIL get_jacobian: onames not an Array:", jac_call_log);
      fail++;
    } else if (jac_call_log.opts_type !== "object") {
      console.log("FAIL get_jacobian: opts not an object:", jac_call_log);
      fail++;
    } else {
      console.log("ok -- get_jacobian dispatch + arg types: " + JSON.stringify(jac_call_log));
      pass++;
    }
  } catch (e) {
    console.log("FAIL get_jacobian: " + e.message);
    fail++;
  }

  // --- 2. Function.stats() returns a native JS object ---
  // Pre-jsout-typemap fix, this returned the raw EM_VAL handle as a
  // number (e.g. 46) -- broken since pre-refactor wasm-js never
  // implemented from_ptr<map<string, M>>.  Now stats is a plain JS
  // dict, with primitive-typed values via GenericType variant dispatch.
  try {
    const x = M.SX.sym("x");
    const nlp = { x: x, f: M.times(M.minus(x, M.SX(2)), M.minus(x, M.SX(2))) };
    const opts = { ipopt: { print_level: 0, sb: "yes" }, print_time: false };
    const solver = M.nlpsol("solver_stats_smoke", "ipopt", nlp, opts);
    solver.call({ x0: M.DM(0) });
    const s = solver.stats();
    if (typeof s !== "object") throw new Error("stats not object: " + typeof s);
    if (Object.keys(s).length === 0) throw new Error("stats empty");
    if (typeof s.return_status !== "string") throw new Error("return_status not string");
    if (typeof s.success !== "boolean") throw new Error("success not boolean");
    if (typeof s.iter_count !== "number") throw new Error("iter_count not number");
    console.log("ok -- Function.stats() returns native JS object with typed primitive values");
    pass++;
  } catch (e) {
    console.log("FAIL Function.stats(): " + e.message);
    fail++;
  }

  // --- 3. Function instances are directly callable + factory-style construction ---
  // M.Function('name', ...) (preferred, no `new`) AND new M.Function(...)
  // both produce the same callable.  `f([args])` is sugar for `f.call([args])`.
  // Methods like f.expand() also return callables (post-wrapped to handle the
  // ES6 inner-class-binding immutability trap -- see casadi.i %insert("js")).
  try {
    const x = M.SX.sym("x");
    // Factory style (preferred, matches existing test convention).
    const f = M.Function("sq", [x], [M.times(x, x)]);
    if (typeof f !== "function") throw new Error("factory: typeof f not function");
    if (!(f instanceof M.Function)) throw new Error("factory: not instanceof M.Function");
    // `new`-style still works for backward compat.
    const f_new = new M.Function("sq_new", [x], [M.times(x, x)]);
    if (typeof f_new !== "function") throw new Error("new: typeof f_new not function");
    // Callable (sugar, n_out-dependent return) -- here n_out=1 so bare DM.
    // Use variadic-positional: one arg per C++ parm.  Wrapping in a JS
    // array `f([M.DM(3)])` would pass the WHOLE array as input 1 (auto-
    // coercing to DM via casadi's array-to-DM rule for primitives, or
    // erroring for a DM-of-DM).
    if (Number(f(M.DM(3)).nonzeros()[0]) !== 9) throw new Error("f(DM(3)) != 9");
    // f.call (programmatic, ALWAYS returns a list).
    const callList = f.call([M.DM(4)]);
    if (!Array.isArray(callList)) throw new Error("f.call must return list, got " + typeof callList);
    if (Number(callList[0].nonzeros()[0]) !== 16) throw new Error("f.call([4])[0] != 16");
    // Variadic with a primitive arg also works (auto-coerce via DM).
    if (Number(f(5).nonzeros()[0]) !== 25) throw new Error("f(5) variadic != 25");
    // Dict-style (input names, not symbol names -- anonymous SX input -> "i0")
    const r3 = f({ i0: M.DM(5) });
    if (typeof r3 !== "object" || Number(r3.o0.nonzeros()[0]) !== 25) throw new Error("dict call != 25");
    // Function returned from a method must ALSO be callable (covered via
    // post-wrap of Function.prototype methods)
    const g = f.expand();
    if (typeof g !== "function") throw new Error("f.expand() not function");
    if (Number(g(M.DM(6)).nonzeros()[0]) !== 36) throw new Error("g(6) != 36");
    console.log("ok -- M.Function(...) factory + new M.Function(...); f([x]) direct call; f.expand() returns callable");
    pass++;
  } catch (e) {
    console.log("FAIL callable Function: " + e.message);
    fail++;
  }

  console.log(`\n${pass}/${pass+fail} passed`);
  process.exit(fail === 0 ? 0 : 1);
})().catch(e => { console.error("FATAL:", e); process.exit(1); });
