# JS Test Suite — Python Port Status

This directory holds a port of selected `test/python/*.py` modules to
JavaScript, targeting the wasm-js build (`build-wasm/swig/wasm-js/casadi.js`).

## Coverage Summary

| Python module          | JS file                | Tests ported | Total Python tests | Status |
|------------------------|------------------------|-------------:|-------------------:|--------|
| `mx.py`                | `mx.js`                | 47           | (43 covers core)   | Substantive |
| `ad.py`                | `ad.py` → `ad.js`      | 9            | (covers SX+MX AD)  | Substantive |
| `matrix.py`            | `matrix.js`            | 18           | (covers core)      | Substantive |
| `multiplication.py`    | `multiplication.js`    | 4            |                    | Substantive |
| `sparsity.py`          | `sparsity.js`          | 18           | 26                 | Substantive |
| `conic.py`             | `conic.js`             | 3            | (highs path)       | Smoke |
| `function.py`          | `function.js`          | 38           | 117                | Substantive |
| `sx.py`                | `sx.js`                | 24           | 80                 | Subset |
| `misc.py`              | `misc.js`              | 9            | 21                 | Subset |
| `vectortools.py`       | `vectortools.js`       | 1            | 1                  | Complete |
| `nlp.py`               | `ipopt_smoke.js` + `ipopt_callback_smoke.js` | 2 | (ipopt only) | Smoke |
| `optistack.py`         | `optistack.js`         | 4            | 20                 | Subset |
| `serialize.py`         | `serialize.js`         | 7            | 3                  | **Substantive** |
| `integration.py`       | `integration.js`       | 6            | 30+ (rk-only)      | Subset |

Director-side tests (no Python equivalent — wasm-js-specific):
| JS file                      | Tests | Purpose                                   |
|------------------------------|------:|-------------------------------------------|
| `director.js`                | 3     | Basic Callback director dispatch          |
| `director_extras_smoke.js`   | 3     | get_jacobian, stats(), callable Function  |
| `phase0_dispatch.js`         | 6     | MX-from-string/null-arg type checks       |

**Total: 13 JS test files, 187 tests passing.**

## Skipped Modules (with rationale)

These Python test modules were not ported because they require infrastructure
not available in the wasm-js build, or depend on Python-specific helpers /
libraries that have no JS equivalent.

| Python module        | Reason for skip                                                                 |
|----------------------|---------------------------------------------------------------------------------|
| `typemaps.py`        | Numpy / scipy.sparse interop -- no JS analog (40 tests)                          |
| `tools.py`           | `casadi.tools.struct` / `struct_symSX` / `entry` (Python-only) (22 tests)        |
| `integration.py`     | Requires CVODES integrator plugin (not in wasm build)                            |
| `simulator.py`       | Requires CVODES integrator (2 tests)                                             |
| `nlp.py`             | Requires Ipopt + other nlpsol plugins; basic Ipopt covered via `ipopt_smoke.js`  |
| `ocp.py`             | Optimal control problems — requires multiple plugins                             |
| `implicitfunction.py`| `rootfinder` API + Newton/Kinsol plugins; not enough wasm-side support           |
| `linearsolver.py`    | `Linsol` plugins (csparse / lapack-lu); not enabled in wasm build                |
| `threads.py`         | Multi-threading (Function::map(..., 'thread', ...))                             |
| `daebuilder.py`      | DaeBuilder is exposed but extensively tied to FMI/file I/O                       |
| `pyright_stubs.py`   | Validates Python pyright stubs — not applicable to JS (TypeScript stubs separate)|
| `feasiblesqpmethod.py`| Specific SQP method solver plugin                                                |
| `serialize.py`       | Function::serialize has known wasm-js OOB bug; smoke-tested in misc.js skip      |
| `optistack.py`       | Opti exposed but `new Opti()` ctor dispatch has SWIG-emit string-marshaling gap   |

## Known Pre-existing wasm-js Gaps Surfaced by the Port

These are documented inline in `function.js` and `optistack.js`:

1. ~~`Function::name()` returns `"emsc"`~~ **FIXED.** Same root cause as
   #5; closed by the constructorHandler patch.

2. **`default_in` option ignored by dict-call** -- dict keys missing from
   the call default to `DM(0)`, not the `default_in` value.

3. **Unknown options silently accepted** at construction and eval time.
   Python casadi raises; wasm-js currently does not.

4. ~~`Function::serialize` round-trip → "memory access out of bounds"~~
   **FIXED.** Was caused by Function::name() being garbage (SWIG-emit
   string-marshaling gap) + Function.deserialize picking the wrong
   overload (first-wins-per-arity in `staticmemberfunctionHandler`).
   Both fixed; round-trip works.  See `serialize.js`.

5. ~~Constructor dispatcher passes JS strings raw to wasm exports~~
   **FIXED.** Patched in `Source/Modules/wasm_js.cxx::constructorHandler`:
   string-typed parms (detected via `tmap:ctype == const char*`) now
   emit the `M.lengthBytesUTF8 + M._malloc + M.stringToUTF8 + M._free`
   dance.  Closes Function::name(), `new Opti(...)`, and many other
   paths that pass string args via the constructor dispatcher.  Also
   fixed a pre-existing related bug: the `prim_checks` typeof-string
   gate used to false-positive on `std::map<std::string, T>` (substring
   match on "string"); now uses ctype-discrimination too.

6. **`Function::size_in(i)` returns `pair<int,int>`** — no JS typemap.
   Workaround: use `size1_in` / `size2_in` / `nnz_in`.

## Fixed During the Port

7. ~~`Function.prototype.call(...)` was DM-only~~. Fixed: now dispatches
   to vector<DM>/<MX>/<SX> + DMDict/SXDict/MXDict based on input type.
   See `casadi.i %insert("js")` and `test_call_overload_dispatch` in
   `function.js`.

8. ~~Constructor dispatcher string-arg marshaling gap~~ FIXED (see #5
   above).

## Remaining Newly-Surfaced Gaps

9. ~~`Function::deserialize(string)` picks the wrong overload~~ **FIXED.**
   Same root cause as the constructor first-wins-per-arity trap; fixed
   by adding type_checks-based discrimination to
   `staticmemberfunctionHandler` AND `memberfunctionHandler` dispatchers
   in `Source/Modules/wasm_js.cxx`.  Member methods like
   `SerializerBase::pack(DM|MX|SX|Sparsity)` now also dispatch on arg
   type instead of silently picking the first-emitted overload.

10. ~~Abstract base classes without their own JS constructor~~ **FIXED.**
    Patched `classHandler` to emit a stub constructor (PRIVATE_CTOR
    only) when no public ctors are collected.  Affects:
    `StringSerializer`, `StringDeserializer`, `SerializerBase`, etc.

11. ~~`Opti.variable()` 0-arg dispatch missing~~ **FIXED.**
    `memberfunctionHandler` now emits default-arg phantom truncated
    overloads (mirrors `staticmemberfunctionHandler`).  Plus also
    fixed:
    - `build_defaults_prologue` BigInt-vs-Number heuristic was too
      eager (matched any "int" substring → produced BigInt for plain
      `int` parms, causing "Cannot convert BigInt to number" on i32
      wasm exports).  Now checks for `long`/`casadi_int`/`int64`
      markers first; plain `int` falls through to plain JS number.
    - `build_defaults_prologue` previously emitted `undefined` for
      C++ ctor-call defaults like `Dict()` / `MXVector()`; the wasm
      side then errored "Failed to convert input N to type '...'".
      Now emits `{}` for Dict-shaped, `[]` for vector-shaped,
      `new <Class>()` for registered class types.
    - `OptiSol::value()` etc. returned `native_DM` which routed
      through `full_or_sparse` -- a helper not implemented for
      wasm-js (returns 0).  Added a SWIGWASMJS-gated typemap
      override that uses `casadi::from_ref` like the regular DM
      out path.

## Running

```sh
cd test/javascript
for f in *.js; do node $f; done
```

Or with the per-build script wrapper used in CI.
