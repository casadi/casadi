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
  const ca = await create();                 // async: instantiates the wasm module
  const x = ca.SX.sym('x', 2);
  const [x0, x1] = ca.vertsplit(x);
  const nlp = { x, f: ca.plus(ca.times(x0, x0), ca.times(x1, x1)),
                   g: ca.minus(ca.plus(x0, x1), ca.SX(10)) };

  await ca.load_nlpsol('ipopt');             // on-demand plugin load (async)
  const solver = ca.nlpsol('solver', 'ipopt', nlp);
  const sol = solver.call({ lbg: ca.DM(0) });
  console.log(sol['x'].nonzeros());         // [5, 5]
})();
```

JavaScript has no operator overloading — use `ca.plus`, `ca.minus`,
`ca.times`, etc. instead of `+ - *`. Solver plugins
(`ipopt`, `fatrop`, sundials integrators, ...) ship as sibling `.so`
side modules and load on demand via `await ca.load_<type>('<name>')`.

See `examples/` for ports of the CasADi documentation examples.
