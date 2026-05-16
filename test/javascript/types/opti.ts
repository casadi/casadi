// opti.ts -- Type-level coverage of the Opti stack (NLP problem builder).

import { Opti, OptiSol, MX, DM, Sparsity, sumsqr, plus, minus, times, _GenericType } from "./casadi";

function expectType<T>(_v: T): void {}

// ============================================================
// Construct an NLP via Opti
// ============================================================
const opti = new Opti();
expectType<Opti>(opti);

// variable/parameter return MX
const x = opti.variable();
expectType<MX>(x);
const y = opti.variable(2n);
expectType<MX>(y);
const z = opti.variable(3n, 3n);
expectType<MX>(z);
const z_sp = opti.variable(Sparsity.lower(3n));
expectType<MX>(z_sp);

const p = opti.parameter();
expectType<MX>(p);
const p2 = opti.parameter(4n, 1n);
expectType<MX>(p2);

// minimize: returns void (sets the cost)
opti.minimize(sumsqr(minus(x, p)));

// subject_to with single constraint and vector of constraints
opti.subject_to(plus(x, p));
opti.subject_to([plus(x, p), minus(y, p2)]);
// Note: Opti.subject_to does NOT have a dict-form overload (unlike
// nlpsol input/output dicts).  Constraints are unnamed positional.
opti.subject_to();  // clears

// solver setup
opti.solver("ipopt");
opti.solver("ipopt", { expand: true });
opti.solver("ipopt", { expand: true }, { print_level: 0, sb: "yes" });

// solve returns an OptiSol
const sol = opti.solve();
expectType<OptiSol>(sol);
expectType<DM>(sol.value(x));
expectType<Record<string, _GenericType>>(sol.stats());

// value() returns DM regardless of input expression type
expectType<DM>(opti.value(x));
expectType<DM>(opti.value(times(x, p)));
expectType<DM>(opti.value(sumsqr(x)));

void [opti, x, y, z, z_sp, p, p2, sol];
