# casadi + pint units (experiment)

Tested with casadi 3.8.1, pint 0.25.3 and numpy 2.4 (`pip install casadi pint numpy`).

| script | what |
|---|---|
| `probe_raw_pint.py` | Experiments 1 and 2: a raw SX/MX, or a `ca.ArrayInterface`, inside a `pint.Quantity` |
| `casadi_pint.py`    | Experiment 3: the `CQ` helper class (casadi expression + pint unit, with `__SX__`/`__MX__`/`__DM__`) |
| `demo_casadi_pint.py` | Runs `CQ` checks and a small optimal-control problem posed in mixed units |

## 1. Raw SX/MX in a pint Quantity: arithmetic works, casadi functions strip units

`ureg.Quantity(ca.SX.sym('x', 3), 'm')` (or `x * ureg.m`) works for SX and MX.
These all work with correct units: `+ - * / **`, unit conversion (`q + 3*ureg.cm`, `.to('km')`),
indexing, `np.sin(q/ureg.m)`, `np.sqrt(q)`, `np.sum(q)`, and `@` with numpy arrays.
Mismatched dimensions raise `DimensionalityError`.

**The catch:** every casadi-namespace function silently drops the unit.
Affected: `ca.sin(q)` with `q` in metres, `ca.vertcat`, `ca.jacobian`, `ca.mtimes` and `ca.dot`.
It happens because pint's `Quantity.__getattr__` forwards SWIG's `.this` pointer
to the magnitude, so SWIG treats the Quantity as the bare SX/MX.
Comparisons (`q <= 5*ureg.m`) return plain casadi expressions after conversion to a common unit.
That is what a constraint needs.
`ca.Function('f', [q], ...)` fails.

## 2. `ca.ArrayInterface` inside pint: same results, plus one extra failure

`Quantity(ArrayInterface(x), 'm')` and `ureg.m * ArrayInterface(x)` work, with the same results as case 1.
`ArrayInterface(x) * ureg.m` raises `DimensionalityError`.
The cause: `ArrayInterface.__mul__` (`__array_priority__` 1005) accepts the `Unit` itself instead of returning `NotImplemented`.
The array interface does not solve the unit-stripping problem.

## 3. `CQ` helper class

- `CQ` holds a `pint.Quantity` whose magnitude is a casadi expression.
  pint does all the unit algebra.
- `CQ` has no `__getattr__` forwarding, so casadi cannot strip the unit without you noticing.
- `__SX__`/`__MX__`/`__DM__` let a `CQ` go straight into casadi (`ca.Function`, `ca.vertcat`, `Opti`, ...).
  The expression is passed in the CQ's own unit; dimensionless ratios are simplified first, so `x[m]/1cm` becomes `100*x`.
  With `CQ.strict = True`, only dimensionless quantities are converted.
  A dimensional `ca.sin(x)` is then rejected, but SWIG replaces the unit error with casadi's generic "wrong type" error.
- Unit-aware helpers: `cp.sin/exp/log/...` (dimensionless argument only), `sqrt`, `fabs`, `sumsqr`, `sum1`, `vertcat/horzcat` (convert to the first argument's unit), `jacobian/gradient` (unit y/x), `mtimes`, and a `Function` wrapper that returns the in/out units.
- Comparisons (`==`, `<=`, ...) convert the right-hand side to the left-hand unit and return a casadi expression, ready for `opti.subject_to`.
  pint's own `__eq__` calls `bool()` on the magnitude, which fails for MX.
- `CQ` is registered as a pint "upcast type", so `Quantity <op> CQ` hands over to CQ.

The demo solves a small optimal-control problem with IPOPT, stated in mixed units:
speed in km/h, acceleration in m/s², time step in minutes, distance in km.
All conversions happen automatically.

Known limits (proof of concept only):
- A CQ carries one unit for the whole matrix.
  A state vector with mixed units needs one CQ per block.
- `CQ` checks units only while the expression is built; the casadi graph itself has no unit information.
