// sparsity.ts -- Type-level coverage of the Sparsity API.

import { Sparsity, DM, SX, MX } from "./casadi";

function expectType<T>(_v: T): void {}

// ============================================================
// Factory methods
// ============================================================
expectType<Sparsity>(Sparsity.dense(3, 4));
expectType<Sparsity>(Sparsity.dense(3));            // square shorthand
expectType<Sparsity>(Sparsity.dense([3, 4]));      // [bigint, bigint] form
expectType<Sparsity>(Sparsity.scalar());
expectType<Sparsity>(Sparsity.scalar(true));
expectType<Sparsity>(Sparsity.unit(5, 2));
expectType<Sparsity>(Sparsity.upper(4));
expectType<Sparsity>(Sparsity.lower(4));
expectType<Sparsity>(Sparsity.diag(3));
expectType<Sparsity>(Sparsity.diag(3, 5));
expectType<Sparsity>(Sparsity.banded(5, 1));
expectType<Sparsity>(Sparsity.rowcol([0, 1], [0, 2], 3, 3));
expectType<Sparsity>(Sparsity.triplet(3, 3, [0, 1], [0, 2]));
// triplet with invert_mapping returns tuple [Sparsity, mapping]
const [sp_tri, mapping] = Sparsity.triplet(3, 3, [0, 1], [0, 2], true);
expectType<Sparsity>(sp_tri);
expectType<bigint[]>(mapping);

// ============================================================
// Introspection
// ============================================================
const sp = Sparsity.dense(3, 4);
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
const sp1 = Sparsity.lower(3);
const sp2 = Sparsity.upper(3);
expectType<Sparsity>(sp1.unite(sp2));
expectType<Sparsity>(sp1.intersect(sp2));
expectType<boolean>(sp1.is_subset(sp2));

// ============================================================
// Matrix-typed triplets
// ============================================================
expectType<DM>(DM.triplet([0, 1], [0, 2], DM([1.0, 2.0])));
expectType<DM>(DM.triplet([0, 1], [0, 2], DM([1.0, 2.0]), 3, 3));
expectType<SX>(SX.triplet([0, 1], [0, 2], SX.sym("v", 2)));

// ============================================================
// Serialization
// ============================================================
expectType<string>(sp.serialize());

void [sp_tri, mapping, sp1, sp2];
