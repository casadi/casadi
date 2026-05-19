# CasADi runtime-helper protocol — canon

This document is the **authoritative contract** for runtime helpers that
live in `casadi/core/runtime/casadi_<name>.hpp` (or in
`casadi/interfaces/<name>/<name>_runtime.hpp` for plugin-private helpers).
It pairs with the diagram in `flow.dot`.

Helpers must be callable from **two paths**:

1. The **C++ vm path**: a class's `init()` / `set_work()` / `eval()`
   chain in casadi_core / casadi_solvers / casadi_interfaces calls the
   helper directly.
2. The **codegen path**: `Function::generate(...)` emits C that calls
   exactly the same helper functions inline.

The protocol exists so the helper is written once and serves both. It
also "avoids codegen/vm duplication of populating helper problem
struct elements" (see `flow.dot`).


## The five canonical functions

A helper named `casadi_<x>` exposes exactly these C runtime functions
(nothing more, nothing less):

| Function | Signature | When called |
|---|---|---|
| `casadi_<x>_setup` | `(prob* p)` | once at C++ class construction (via `set_<x>_prob`); also emitted once per call in codegen |
| `casadi_<x>_work` | `(const prob* p, sz_arg, sz_res, sz_iw, sz_w)` | once at C++ class init (in `alloc_*` block) |
| `casadi_<x>_set_work` | `(data* d, arg, res, iw, w)` | per call; mirrors `FunctionInternal::set_work` |
| `casadi_<x>_solve` (or `_eval`) | `(data* d, ...)` | per call; the numerical kernel |
| `casadi_<x>_init` *(optional)* | `(data* d, iw, w)` | only if the helper has a scratch-only sub-helper used **internally** by other helpers; never emitted from codegen |

### Strict rules

1. **`set_work` signature is `(d, arg, res, iw, w)`.** No extra ints,
   flags, or sizes.  All immutable per-instance configuration goes on
   `prob` (filled in `set_<x>_prob` from class members, or computed in
   `setup`).
2. **`work` signature is `(p, sz_arg, sz_res, sz_iw, sz_w)`.** Same
   reason — anything else needed lives on `prob`.
3. **`set_work` claims ALL per-call buffers from `*iw` / `*w`** and
   assigns to `d->*` pointer fields.  After this call the caller uses
   `d->RSQ_val`, `d->h_hat_csc`, etc. exclusively.  **No raw `*w += N`
   in the codegen body or in the C++ `set_work` override.**
4. **Codegen emits `_set_work` only.** `casadi_<x>_init` is reserved
   for nested internal use (e.g. one helper's `set_work` calling another
   helper's `init` for an embedded scratch claim).  Examples:
   `casadi_socp_init`, `casadi_ipqp_init`, `casadi_jac_init`.
5. **Plugin-private runtime helpers** (in `casadi/interfaces/<name>/`)
   follow the same rules.  See e.g. `casadi_daqp_set_work`,
   `casadi_xpress_set_work`.


## Lifecycle

```
init  (C++ class, called once at Function construction)
 ├─ set_<x>_prob()              -- fill prob from class-immutable fields
 │    └─ casadi_<x>_setup(p)    -- compute derived prob fields
 └─ alloc_w/alloc_iw via casadi_<x>_work(p, &sz_iw, &sz_w)

call/eval  (per call)
 ├─ setup()
 │    └─ set_work(...)
 │        └─ casadi_<x>_set_work(d, &arg, &res, &iw, &w)
 │              i.e. d->* = *w; *w += ...   (claim from scratch)
 └─ solve(...)  / casadi_<x>_solve(d, ...)
```

The codegen emission walks the same shape: it emits a function-local
`prob` struct, calls `casadi_<x>_setup` on it, claims local `data`,
calls `casadi_<x>_set_work`, then `casadi_<x>_solve`.


## What goes on `prob` vs `data`

The split, in one line: **scalars on `prob`, `casadi_int*` arrays on `data`.**

- **`prob`**: read-only after class init. Stores:
  - Caller-owned descriptors (`const`-pointer fields like `sp_a`,
    sparsity pointers, the original partition arrays).
  - **Scalar** setup-derived fields filled by `casadi_<x>_setup`
    (cumulative sizes, max dims, nnz totals).
  - **Solver configuration** (immutable per-instance ints/flags).
    These used to be stuffed onto `set_work` / `work` signatures as
    extra parameters; that's now banned — see `casadi_sqpmethod_prob`'s
    `elastic_mode` / `so_corr` and `casadi_feasiblesqpmethod_prob`'s
    `sz_anderson_memory` for the migration pattern.

- **`data`**: per-call workspace. Pointer fields (including all
  `casadi_int*` derived per-stage arrays) wired by
  `casadi_<x>_set_work` from `*iw`/`*w`. Plus a back-reference to prob.

**Hard rule on derived `casadi_int*` arrays:** if a setup function
needs to write into a `casadi_int*` field at runtime (e.g. per-stage
`nx_hat[]`, `*_offsets[]`, struct-of-block arrays), that field
**belongs on `data`**, not `prob`.  `set_work` claims the iw slots and
populates them via the same per-stage recurrence.  Scalars (`nx_max`,
`nnz_*`, `total_*`) stay on `prob`.

This is enforced because:
- `prob` is conceptually "what the caller passes in plus what setup
  computes scalarwise."  Caller-owned arrays stay caller-owned (const);
  derived arrays would otherwise need *caller-allocated, helper-filled*
  storage on `prob`, which leaks workspace concerns into the prob
  contract and breaks codegen (no class-member backing).
- `data` is "per-call scratch the kernel writes to."  Per-stage
  derived `casadi_int*` arrays are exactly that — they need backing
  storage that lives across set_work → solve → reset.

### Prior art for arrays-on-data

- **`casadi_hpipm_set_work`** populates `d->hlbx[]`, `d->hubx[]`,
  `d->hidxbx[]`, `d->hidxbu[]`, etc. — per-stage pointer/index arrays
  derived from `p->nx[]`, `p->nu[]`, `p->nbx[]`.  Filled in 5+ loops
  after the iw/w claims.
- **`casadi_condensing_set_work`** populates `d->nx_hat[]`,
  `d->nu_hat[]`, `d->ng_hat[]`, `d->AB_hat[]`, `d->CD_hat[]`,
  `d->RSQ_hat[]`, `d->{AB,CD,RSQ}_hat_offsets[]` — all per-condensed-
  stage derived from `p->nx/nu/ng/M` via the condensing recurrence.

In both cases `setup(p)` computes only scalar metadata; `set_work(d, …)`
does the per-stage work that materializes the `casadi_int*` arrays.

### Workspace vs BSS — what stays in BSS

A small amount of prob-side data still legitimately lives in BSS
(function-scope `static`) in codegen output: fixed-size block-descriptor
arrays that are an *input* to the algorithm (not a derivation), filled
once via `casadi_unpack_ocp_blocks` from a length-prefixed
`g.constant()` const.  The condensing codegen body still uses this
pattern for the original (non-hat) `cond_AB_blk`, `cond_CD_blk`,
`cond_RSQ_blk`:

```c
static struct casadi_ocp_block cond_AB_blk[3];
casadi_unpack_ocp_blocks(cond_AB_blk, casadi_s7);
p_cond.AB = cond_AB_blk;
```

These hold the *original problem* descriptors — caller config, not
derived per-call.  The hat counterparts (`AB_hat[]`, etc.), which used
to also be in BSS, now live on data via `iw` since they are derived per
call by the recurrence in `set_work`.

**The line, refined:** prob-side data is BSS-OK only when (a) it's a
direct unpack of a compile-time const and (b) the algorithm consumes it
as a fixed input.  Anything `setup` would otherwise *derive* into a
`casadi_int*` array goes on `data` via `iw`, populated in `set_work`.


## The trivial-vs-non-trivial size rule

When `casadi_<x>_set_work` needs a **size** to claim from `*w`, where
does it get it from?

- **Trivially derivable**: O(1) — direct field access or arithmetic on
  prob fields with no loops. Compute inline.
  `casadi_int n = p->nx[0] + p->nu[0];`
- **Non-trivially derived**: anything with a loop, a sparsity walk, or
  partition logic. Cache it on `prob` (as a *scalar*) and fill it in
  `casadi_<x>_setup`.  `set_work` reads `p->nnz_RSQ` directly.

Even a single-loop sum like `sum_k (nx[k]+nu[k])^2` is non-trivial —
the resulting *scalar* goes on `prob`. The bar is O(1).

Note this rule is about *sizes/scalars*. `casadi_int*` arrays follow
the rule above (arrays-on-data, populated in `set_work`).


## Example references in the casadi tree

- **`casadi_qp`** (`casadi/core/runtime/casadi_qp.hpp`): minimal example
  — `set_work` is empty (QP just hands pointers through).
- **`casadi_condensing`** (`casadi/core/runtime/casadi_condensing.hpp`):
  full example — `setup` computes ~17 cumulative-size scalars onto
  prob; `set_work` claims iw slots for the 9 per-condensed-stage
  derived arrays (`nx_hat`, `nu_hat`, `ng_hat`, `AB_hat`, `CD_hat`,
  `RSQ_hat`, `*_hat_offsets`) AND populates them via the per-K
  recurrence; then claims ~25 w buffers.  The Conic-base codegen body
  collapses to <30 lines of `d_cond.<field>` use.
- **`casadi_hpipm`** (`casadi/interfaces/hpipm/hpipm_runtime.hpp`):
  precedent for the per-stage-arrays-on-data pattern — `set_work`
  populates `d->hlbx[]`, `d->hubx[]`, `d->hidxbx[]`, `d->hidxbu[]`,
  `d->pis[]` via 5+ post-claim loops over `p->nx[]`, `p->nu[]`,
  `p->nbx[]`, `p->nbu[]`.
- **`casadi_oracle`** / **`casadi_nlpsol`**: `set_work` does both
  scratch claim AND arg/res snapshot (oracle) or arg/res bump
  (nlpsol) for sub-Function slot allocation.
- **`casadi_socp`** / **`casadi_ipqp`** / **`casadi_jac`**: examples of
  the `_init` (3-arg) pattern — internal scratch-only helpers, not
  emitted from codegen.


