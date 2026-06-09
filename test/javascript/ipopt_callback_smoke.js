// Test iteration_callback under Node to isolate browser from C++ issues.
const path = require("path");
const create = require(path.resolve(
  "/home/jgillis/programs/casadi/build-wasm-flang/build/swig/wasm-js/casadi.js"
));

(async () => {
  const M = await create();
  console.log("has ipopt:", M.has_nlpsol("ipopt"));

  const x = M.SX.sym("x", 1, 1);
  const y = M.SX.sym("y", 1, 1);
  const v = M.vcat([x, y]);
  const f = M.plus(
    M.times(M.minus(x, M.SX(1)), M.minus(x, M.SX(1))),
    M.times(M.minus(y, M.SX(2)), M.minus(y, M.SX(2)))
  );
  const g = M.plus(x, y);

  const n = 2, m = 1;
  class IterCB extends M.Callback {
    constructor() {
      super();
      this.iter = 0;
      this.construct("IterCB", null);
    }
    get_n_in()  { return 6n; }
    get_n_out() { return 1n; }
    get_name_in(i)  { return ["x","f","g","lam_x","lam_g","lam_p"][Number(i)]; }
    get_name_out(i) { return "ret"; }
    get_sparsity_in(i) {
      const j = Number(i);
      if (j === 0) return M.Sparsity.dense(BigInt(n));
      if (j === 1) return M.Sparsity.dense(1n);
      if (j === 2) return M.Sparsity.dense(BigInt(m));
      if (j === 3) return M.Sparsity.dense(BigInt(n));
      if (j === 4) return M.Sparsity.dense(BigInt(m));
      if (j === 5) return M.Sparsity(0n, 0n);
      return M.Sparsity(0n, 0n);
    }
    get_sparsity_out(i) { return M.Sparsity.dense(1); }
    eval(args) {
      const xs = args[0];
      console.log("  iter", this.iter, "x=", xs[0].nonzeros(), "f=", xs[1].nonzeros());
      this.iter++;
      return [M.DM(0)];
    }
  }
  const cb = new IterCB();

  const opts = {
    ipopt: { print_level: 0, sb: "yes" },
    print_time: false,
    iteration_callback: cb,
  };
  console.log("constructing solver...");
  try {
    const solver = M.nlpsol("solver", "ipopt", { x: v, f: f, g: g }, opts);
    console.log("solving...");
    const res = solver.call({
      x0: M.DM([0, 0]),
      lbx: M.DM([-1e20, -1e20]),
      ubx: M.DM([1e20, 1e20]),
      lbg: M.DM([1]),
      ubg: M.DM([1e20]),
    });
    console.log("x*=", res.x.nonzeros(), "f*=", res.f.nonzeros());
  } catch (e) {
    console.error("CAUGHT:", e);
    if (e instanceof WebAssembly.Exception) {
      // Try to extract the C++ exception message
      const M2 = await create();
      const tags = e.getArg ? "" : "(no getArg method)";
      console.error("tags:", tags);
    }
  }
})();
