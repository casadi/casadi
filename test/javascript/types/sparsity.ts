// sparsity.ts -- Type-level coverage of the Sparsity API.

import { Sparsity, DM, SX, MX } from "./casadi";

function expectType<T>(_v: T): void {}

// ============================================================
// Factory methods
// ============================================================
expectType<Sparsity>(Sparsity.dense(3n, 4n));
expectType<Sparsity>(Sparsity.dense(3n));            // square shorthand
expectType<Sparsity>(Sparsity.dense([3n, 4n]));      // [bigint, bigint] form
expectType<Sparsity>(Sparsity.scalar());
expectType<Sparsity>(Sparsity.scalar(true));
expectType<Sparsity>(Sparsity.unit(5n, 2n));
expectType<Sparsity>(Sparsity.upper(4n));
expectType<Sparsity>(Sparsity.lower(4n));
expectType<Sparsity>(Sparsity.diag(3n));
expectType<Sparsity>(Sparsity.diag(3n, 5n));
expectType<Sparsity>(Sparsity.banded(5n, 1n));
expectType<Sparsity>(Sparsity.rowcol([0n, 1n], [0n, 2n], 3n, 3n));
expectType<Sparsity>(Sparsity.triplet(3n, 3n, [0n, 1n], [0n, 2n]));
// triplet with invert_mapping returns tuple [Sparsity, mapping]
const [sp_tri, mapping] = Sparsity.triplet(3n, 3n, [0n, 1n], [0n, 2n], true);
expectType<Sparsity>(sp_tri);
expectType<bigint[]>(mapping);

// ============================================================
// Introspection
// ============================================================
const sp = Sparsity.dense(3n, 4n);
expectType<bigint>(sp.size1());
expectType<bigint>(sp.size2());
expectType<bigint>(sp.numel());
expectType<bigint>(sp.nnz());
expectType<[bigint, bigint]>(sp.size());
expectType<boolean>(sp.is_dense());
expectType<boolean>(sp.is_scalar());
expectType<boolean>(sp.is_vector());
expectType<boolean>(sp.is_square());
expectType<boolean>(sp.is_symmetric());

// ============================================================
// Set ops
// ============================================================
const sp1 = Sparsity.lower(3n);
const sp2 = Sparsity.upper(3n);
expectType<Sparsity>(sp1.unite(sp2));
expectType<Sparsity>(sp1.intersect(sp2));
expectType<boolean>(sp1.is_subset(sp2));

// ============================================================
// Matrix-typed triplets
// ============================================================
expectType<DM>(DM.triplet([0n, 1n], [0n, 2n], new DM([1.0, 2.0])));
expectType<DM>(DM.triplet([0n, 1n], [0n, 2n], new DM([1.0, 2.0]), 3n, 3n));
expectType<SX>(SX.triplet([0n, 1n], [0n, 2n], SX.sym("v", 2n)));

// ============================================================
// Serialization
// ============================================================
expectType<string>(sp.serialize());

void [sp_tri, mapping, sp1, sp2];
