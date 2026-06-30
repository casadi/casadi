// Quick smoke test: solve a tiny NLP using IPOPT via the wasm build.
//   minimize  (x-1)^2 + (y-2)^2
//   subject to x + y >= 1
// Optimum: x=1, y=2, f=0.
const path = require("path");
const create = require(path.resolve(
  "/home/jgillis/programs/casadi/build-wasm-flang/build/swig/wasm-js/casadi.js"
));

(async () => {
  const M = await create();
  console.log("Has Ipopt plugin:", M.has_nlpsol ? M.has_nlpsol("ipopt") : "(no has_nlpsol)");

  // Use SX (the wasm-js dispatcher picks the SX overload of nlpsol).
  const x = M.SX.sym("x", 1, 1);
  const y = M.SX.sym("y", 1, 1);
  const v = M.vcat([x, y]);
  const f = M.plus(
    M.times(M.minus(x, M.SX(1)), M.minus(x, M.SX(1))),
    M.times(M.minus(y, M.SX(2)), M.minus(y, M.SX(2)))
  );
  const g = M.plus(x, y);
  const nlp = { x: v, f: f, g: g };
  const opts = { ipopt: { print_level: 0, sb: "yes" }, print_time: false };
  console.log("Constructing solver...");
  const solver = M.nlpsol("solver", "ipopt", nlp, opts);
  console.log("Solving...");
  const res = solver.call({
    x0:  M.DM([0, 0]),
    lbx: M.DM([-1e20, -1e20]),
    ubx: M.DM([ 1e20,  1e20]),
    lbg: M.DM([1]),
    ubg: M.DM([1e20]),
  });
  console.log("x* =", res.x.nonzeros());
  console.log("f* =", res.f.nonzeros());
})();
