// Required plugin loads and numerical checks shared by Node and browser tests.
async function testWasmPlugins(ca, manifest) {
  const near = (actual, expected, tolerance = 1e-6) => {
    if (!Number.isFinite(actual) || Math.abs(actual - expected) > tolerance) {
      throw Error(`Expected ${expected}, got ${actual}`);
    }
  };
  const tested = [];
  // Test the new solver dependencies before other plugins can supply their symbols.
  const priority = p => p.name === 'fatrop' ? 0 : p.name === 'mumps' ? 1 : 2;
  for (const plugin of [...manifest.plugins].sort((a, b) => priority(a) - priority(b))) {
    if (plugin.excluded) continue;
    const {kind, name} = plugin;
    const loader = kind === 'xmlfile' ? ca.XmlFile.load_plugin : ca['load_' + kind];
    if (typeof loader !== 'function') throw Error(`No loader for ${kind}:${name}`);
    console.log("Testing " + kind + ":" + name);
    await loader(name);
    tested.push(kind + ':' + name);

    if (kind === 'blas') {
      const values = ca.mtimes(ca.DM([[1, 2], [3, 4]]), ca.DM([5, 6]), name).nonzeros();
      near(values[0], 17); near(values[1], 39);
    }
    if (kind === 'xmlfile') ca.XmlFile(name);
    if (kind === 'nlpsol' && ['fatrop', 'ipopt'].includes(name)) {
      const x = ca.SX.sym('x', 2);
      const [a, b] = ca.vertsplit(x);
      const solver = ca.nlpsol(name + '_smoke', name,
        {x, f: ca.dot(x, x), g: ca.plus(a, b)}, name === 'fatrop' ? {print_time: false, fatrop: {tol: 1e-7}}
          : {print_time: false, ipopt: {print_level: 0, sb: 'yes'}});
      const sol = solver.call({x0: ca.DM([0, 0]), lbg: ca.DM(10)});
      if (!solver.stats().success) throw Error(name + " NLP did not converge");
      for (const v of sol.x.nonzeros()) near(v, 5);
      near(sol.f.nonzeros()[0], 50, 1e-4);
    }
    if (kind === 'conic' && ['fatrop', 'osqp'].includes(name)) {
      const x = ca.SX.sym('x', 2);
      const solver = ca.qpsol('qp_' + name, name,
        {x, f: ca.dot(x, x)}, name === 'osqp' ? {osqp: {verbose: false}} :
          {});
      const sol = solver.call({lbx: ca.DM([1, 2]), ubx: ca.DM([10, 10])});
      const values = sol.x.nonzeros();
      near(values[0], 1, 5e-3); near(values[1], 2, 5e-3);
    }
    if (kind === 'rootfinder' && name === 'bisection') {
      const x = ca.SX.sym('root_x');
      const f = ca.Function('residual', [x], [ca.minus(ca.times(x, x), ca.SX(2))]);
      const solver = ca.rootfinder('bisect', name, f, {lb: 0, ub: 2});
      near(solver(ca.DM(1)).nonzeros()[0], Math.sqrt(2));
    }
    if (kind === 'linsol' && ['csparsecholesky', 'mumps'].includes(name)) {
      const a = ca.DM([[4, 1], [1, 3]]);
      const solver = ca.Linsol('linear_' + name, name, a.sparsity());
      const x = solver.solve(a, ca.DM([1, 2])).nonzeros();
      near(x[0], 1 / 11); near(x[1], 7 / 11);
    }
  }
  return tested;
}
if (typeof module !== 'undefined') module.exports = {testWasmPlugins};
