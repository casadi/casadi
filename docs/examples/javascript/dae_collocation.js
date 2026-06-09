//
//     MIT No Attribution
//
//     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//     Permission is hereby granted, free of charge, to any person obtaining a copy of this
//     software and associated documentation files (the "Software"), to deal in the Software
//     without restriction, including without limitation the rights to use, copy, modify,
//     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//     permit persons to whom the Software is furnished to do so.
//
//     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
// JS port of docs/examples/python/dae_collocation.py.
//
// Hand-built direct collocation of a crane/pendulum DAE optimal control
// problem (Mario Zanon & Sebastien Gross, KU Leuven 2012), solved as one
// big NLP with ipopt.  No integrator plugin is used -- the collocation
// equations are assembled by hand -- so this ports directly.  The original
// ipopt option linear_solver='ma27' is dropped (default mumps is used) and
// ipopt.max_iter is left at its default (see the note near the solve).
//
// JS notes (see README.md):
//   * numpy linear algebra (inv, @) is replaced by casadi DM ops (ca.inv,
//     ca.mtimes).  Collocation coefficients C/D are plain JS arrays.
//   * The NLP vector V is sliced via per-scalar vertsplit + range vertcat.
//   * Plotting is dropped; the optimal cost and a few state samples logged.

async function example(ca, log) {
  const inf = Infinity;

  // -------- Collocation setup --------
  const nicp = 1;
  const xref = 0.1;
  const l = 1.0, mmass = 1.0, Mmass = 1.0, g = 9.81;
  const tf = 5.0;
  const nk = 50;
  const ndstate = 6, nastate = 1, ninput = 1;
  const deg = 4;
  const h = tf / nk / nicp;

  const tau = ca.SX.sym("tau");
  const tau_root = [0].concat(ca.collocation_points(deg, "radau"));

  // Lagrange-polynomial collocation coefficients
  const C = Array.from({ length: deg + 1 }, () => new Array(deg + 1).fill(0));
  const D = new Array(deg + 1).fill(0);
  for (let j = 0; j <= deg; ++j) {
    let L = ca.SX(1);
    for (let j2 = 0; j2 <= deg; ++j2) {
      if (j2 !== j) L = ca.times(L, ca.rdivide(ca.minus(tau, tau_root[j2]), tau_root[j] - tau_root[j2]));
    }
    const lfcn = ca.Function("lfcn", [tau], [L]);
    D[j] = Number(lfcn(1.0));
    const tfcn = ca.Function("tfcn", [tau], [ca.tangent(L, tau)]);
    for (let j2 = 0; j2 <= deg; ++j2) C[j][j2] = Number(tfcn(tau_root[j2]));
  }

  // -------- Model setup (implicit DAE) --------
  const t = ca.SX.sym("t");
  const u = ca.SX.sym("u");
  const xd = ca.SX.sym("xd", ndstate);
  const xa = ca.SX.sym("xa", nastate);
  const xddot = ca.SX.sym("xdot", ndstate);
  const p = ca.SX.sym("p", 0, 1);

  const [x, y, w, dx, dy, dw] = ca.vertsplit(xd);
  const [xa0] = ca.vertsplit(xa);
  const xdd = ca.vertsplit(xddot);

  const res = ca.vertcat(
    ca.minus(xdd[0], dx),
    ca.minus(xdd[1], dy),
    ca.minus(xdd[2], dw),
    ca.plus(ca.times(mmass, xdd[3]), ca.times(ca.minus(x, w), xa0)),
    ca.minus(ca.plus(ca.times(mmass, xdd[4]), ca.times(y, xa0)), g * mmass),
    ca.plus(ca.plus(ca.times(Mmass, xdd[5]), ca.times(ca.minus(w, x), xa0)), u),
    ca.plus(ca.plus(ca.plus(
      ca.times(ca.minus(x, w), ca.minus(xdd[3], xdd[5])),
      ca.times(y, xdd[4])),
      ca.times(dy, dy)),
      ca.times(ca.minus(dx, dw), ca.minus(dx, dw))));

  const ffcn = ca.Function("ffcn", [t, xddot, xd, xa, u, p], [res]);
  const MayerTerm = ca.Function("mayer", [t, xd, xa, u, p],
    [ca.plus(ca.plus(ca.plus(
      ca.times(ca.minus(x, xref), ca.minus(x, xref)),
      ca.times(ca.minus(w, xref), ca.minus(w, xref))),
      ca.times(dx, dx)), ca.times(dy, dy))]);
  const LagrangeTerm = ca.Function("lagrange", [t, xd, xa, u, p],
    [ca.plus(
      ca.times(ca.minus(x, xref), ca.minus(x, xref)),
      ca.times(ca.minus(w, xref), ca.minus(w, xref)))]);

  // Bounds
  const u_min = [-2], u_max = [2];
  const xD_min = [-inf, -inf, -inf, -inf, -inf, -inf];
  const xD_max = [inf, inf, inf, inf, inf, inf];
  const xDi_min = [0.0, l, 0.0, 0.0, 0.0, 0.0];
  const xDi_max = [0.0, l, 0.0, 0.0, 0.0, 0.0];
  const xD_init = [0.0, l, 0.0, 0.0, 0.0, 0.0];
  const xA_min = [-inf], xA_max = [inf];
  const xA_init = [Math.sign(l) * 9.81];

  // -------- NLP variable layout --------
  const nx = ndstate + nastate;
  const ndiff = ndstate, nalg = nastate, nu = ninput, NP = 0;
  const NXD = nicp * nk * (deg + 1) * ndiff;
  const NXA = nicp * nk * deg * nalg;
  const NU = nk * nu;
  const NXF = ndiff;
  const NV = NXD + NXA + NU + NXF + NP;

  const V = ca.MX.sym("V", NV);
  const Vrows = ca.vertsplit(V, Array.from({ length: NV + 1 }, (_, i) => i));
  const slice = (a, b) => ca.vcat(Vrows.slice(a, b));

  const vars_lb = new Array(NV).fill(0);
  const vars_ub = new Array(NV).fill(0);
  const vars_init = new Array(NV).fill(0);
  const setRange = (arr, off, vals) => { for (let i = 0; i < vals.length; ++i) arr[off + i] = vals[i]; };

  // XD[k][i][j], XA[k][i][j-1], U[k]
  const XD = Array.from({ length: nk + 1 }, () => Array.from({ length: nicp }, () => new Array(deg + 1).fill(null)));
  const XA = Array.from({ length: nk }, () => Array.from({ length: nicp }, () => new Array(deg).fill(null)));
  const U = new Array(nk).fill(null);

  let offset = 0;
  for (let k = 0; k < nk; ++k) {
    for (let i = 0; i < nicp; ++i) {
      for (let j = 0; j <= deg; ++j) {
        XD[k][i][j] = slice(offset, offset + ndiff);
        if (j !== 0) XA[k][i][j - 1] = slice(offset + ndiff, offset + ndiff + nalg);
        if (k === 0 && j === 0 && i === 0) {
          setRange(vars_init, offset, xD_init);
          setRange(vars_lb, offset, xDi_min);
          setRange(vars_ub, offset, xDi_max);
          offset += ndiff;
        } else if (j !== 0) {
          setRange(vars_init, offset, xD_init.concat(xA_init));
          setRange(vars_lb, offset, xD_min.concat(xA_min));
          setRange(vars_ub, offset, xD_max.concat(xA_max));
          offset += nx;
        } else {
          setRange(vars_init, offset, xD_init);
          setRange(vars_lb, offset, xD_min);
          setRange(vars_ub, offset, xD_max);
          offset += ndiff;
        }
      }
    }
    U[k] = slice(offset, offset + nu);
    setRange(vars_lb, offset, u_min);
    setRange(vars_ub, offset, u_max);
    setRange(vars_init, offset, [0.0]);
    offset += nu;
  }
  XD[nk][0][0] = slice(offset, offset + ndiff);
  setRange(vars_lb, offset, xD_min);
  setRange(vars_ub, offset, xD_max);
  setRange(vars_init, offset, xD_init);
  offset += ndiff;
  if (offset !== NV) throw new Error("offset != NV (" + offset + " vs " + NV + ")");

  const P = ca.MX(0, 1);

  // -------- Constraints --------
  const gcon = [];
  // Collocation + continuity equations
  for (let k = 0; k < nk; ++k) {
    for (let i = 0; i < nicp; ++i) {
      for (let j = 1; j <= deg; ++j) {
        let xp_jk = ca.MX.zeros(ndiff);
        for (let j2 = 0; j2 <= deg; ++j2) xp_jk = ca.plus(xp_jk, ca.times(C[j2][j], XD[k][i][j2]));
        const fk = ffcn(ca.DM(0), ca.rdivide(xp_jk, h), XD[k][i][j], XA[k][i][j - 1], U[k], P);
        const fkrows = ca.vertsplit(fk, [0, ndiff, ndiff + nalg]);
        gcon.push(fkrows[0]); // differential part == 0
        gcon.push(fkrows[1]); // algebraic part == 0
      }
      let xf_k = ca.MX.zeros(ndiff);
      for (let j = 0; j <= deg; ++j) xf_k = ca.plus(xf_k, ca.times(D[j], XD[k][i][j]));
      if (i === nicp - 1) gcon.push(ca.minus(XD[k + 1][0][0], xf_k));
      else gcon.push(ca.minus(XD[k][i + 1][0], xf_k));
    }
  }

  // -------- Objective --------
  // Mayer term at the final collocation node
  let Obj = MayerTerm(ca.DM(0), XD[nk - 1][0][deg], XA[nk - 1][0][deg - 1], U[nk - 1], P);

  // Lagrange term via the collocation quadrature weights
  // lDotAtTauRoot = C.T ; ldInv = inv(C.T[1:,1:]) ; lAtOne = D
  const CT = Array.from({ length: deg + 1 }, (_, a) => Array.from({ length: deg + 1 }, (_, b) => C[b][a]));
  const ldInv = ca.MX(ca.inv(ca.DM(CT.slice(1).map((row) => row.slice(1)))));
  const lAtOne1 = ca.MX(ca.DM(D.slice(1)));

  for (let k = 0; k < nk; ++k) {
    for (let i = 0; i < nicp; ++i) {
      const dqParts = [];
      for (let j = 1; j <= deg; ++j) {
        dqParts.push(LagrangeTerm(ca.DM(0), XD[k][i][j], XA[k][i][j - 1], U[k], P));
      }
      const dQs = ca.times(h, ca.vcat(dqParts)); // deg x 1
      const Qs = ca.mtimes(ldInv, dQs);                       // deg x 1
      Obj = ca.plus(Obj, ca.mtimes(ca.transpose(Qs), lAtOne1)); // scalar
    }
  }
  // -------- Solve --------
  // NB: integer-valued ipopt sub-options (e.g. ipopt.max_iter) are mis-cast
  // by the wasm option bridge and clamp the solve to 1 iteration, so we rely
  // on the default max_iter here.
  const nlp = { x: V, f: Obj, g: ca.vcat(gcon) };
  const solver = ca.nlpsol("solver", "ipopt", nlp,
    { expand: true, "ipopt.tol": 1e-4, "ipopt.print_level": 0, print_time: false });

  const sol = solver.call({
    x0: ca.DM(vars_init), lbx: ca.DM(vars_lb), ubx: ca.DM(vars_ub), lbg: ca.DM(0), ubg: ca.DM(0),
  });

  log("optimal cost: " + sol["f"]);

  // Chariot position x (xd[0]) at the start of each finite element: the
  // per-element layout has a fixed stride, so a strided slice does it.
  const stride = nicp * ((deg + 1) * ndiff + deg * nalg) + nu;
  const xstarts = sol["x"]["0::" + stride];  // nk element starts + final state
  log("chariot position x at element starts (every 10th):");
  log("  " + xstarts["0::10"]);
  log("final x = " + xstarts[nk] + " (target xref = " + xref + ")");
}

if (typeof require !== "undefined" && typeof module !== "undefined" && require.main === module) {
  const path = require("path");
  const casadiPath = process.env.CASADI_JS
    || path.resolve(__dirname, "../../../build-wasm/swig/wasm-js/casadi.js");
  require(casadiPath)()
    .then((ca) => example(ca, (...a) => console.log(...a)))
    .catch((e) => { console.error("FATAL:", e.message || e); process.exit(1); });
}

if (typeof module !== "undefined" && module.exports) module.exports = example;
