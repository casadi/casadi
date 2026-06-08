# @casadi/casadi-wasm

[CasADi](https://web.casadi.org) compiled to WebAssembly. Symbolic
framework for algorithmic differentiation and numerical optimization.

## Install

```
npm install @casadi/casadi-wasm
```

## Use (Node.js)

```js
const create = require('@casadi/casadi-wasm');

(async () => {
  const M = await create();                 // async: instantiates the wasm module
  const x = M.SX.sym('x', 2);
  const [x0, x1] = M.vertsplit(x);
  const nlp = { x, f: M.plus(M.times(x0, x0), M.times(x1, x1)),
                   g: M.minus(M.plus(x0, x1), M.SX(10)) };

  await M.load_nlpsol('ipopt');             // on-demand plugin load (async)
  const solver = M.nlpsol('solver', 'ipopt', nlp);
  const sol = solver.call({ lbg: M.DM(0) });
  console.log(sol['x'].nonzeros());         // [5, 5]
})();
```

JavaScript has no operator overloading — use `M.plus`, `M.minus`,
`M.times`, etc. instead of `+ - *`. Solver plugins
(`ipopt`, `fatrop`, sundials integrators, ...) ship as sibling `.so`
side modules and load on demand via `await M.load_<type>('<name>')`.

See `examples/` for ports of the CasADi documentation examples.
