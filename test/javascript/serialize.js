// serialize.js -- port of test/python/serialize.py (2/3 tests).
//
// Status of casadi's serialize/deserialize machinery on wasm-js:
//
//   WORKING (after the SWIG-emit fixes in this branch):
//     - Function::serialize() / Function.deserialize(string) round-trip.
//     - StringSerializer.pack(...) / .encode().
//     - StringDeserializer.decode() / .unpack_<typed>() per-type.
//
//   STILL BROKEN:
//     - test_compat -- needs `serialize_3.5.5` directory of pre-baked
//       .casadi files (filesystem-dependent).  Not portable.
//     - test_pickling_context -- Python pickle module (not applicable).
//
// Skipped tests (1 of 3):
//   test_compat            -- filesystem-dependent .casadi corpus
//   test_pickling_context  -- Python pickle (not applicable)

const { assertEqual, assertTrue, assertAlmost, checkArray,
        assertArrayAlmostEqual, runTests } = require("./_helpers");

function test_function_serialize_roundtrip(M) {
  // Function::serialize() -> string -> Function::deserialize() preserves
  // name + numeric behaviour.  Bug 1 (Function::name returning "emsc")
  // and Bug 2 (deserialize picking the wrong overload) were both wasm-js
  // SWIG-emit gaps; fixed in Source/Modules/wasm_js.cxx.
  const x = M.SX.sym("x");
  const f = M.Function("squared", [x], [M.times(x, x)]);
  const blob = f.serialize();
  assertTrue(typeof blob === "string" && blob.length > 100,
             "serialize -> non-trivial string");
  const g = M.Function.deserialize(blob);
  assertEqual(g.name(), "squared", "deserialized name preserved");
  // Round-tripped numeric eval matches.
  checkArray(g(3),  [9],  10, "g(3) == 9");
  checkArray(g(-4), [16], 10, "g(-4) == 16");
  // n_in / n_out preserved.
  assertEqual(g.n_in(),  1n, "n_in");
  assertEqual(g.n_out(), 1n, "n_out");
}

function test_string_serialize_dm(M) {
  // StringSerializer / StringDeserializer round-trip for a DM.
  const ss = new M.StringSerializer();
  ss.pack(M.DM([1, 2, 3, 4]));
  const blob = ss.encode();
  assertTrue(typeof blob === "string" && blob.length > 20,
             "StringSerializer blob non-trivial: " + blob.length);
  const sd = new M.StringDeserializer(blob);
  const back = sd.unpack_dm();
  assertTrue(back && back.nonzeros, "unpack_dm() returned DM");
  assertArrayAlmostEqual(back.nonzeros(), [1, 2, 3, 4], 10,
                          "DM round-trip values");
}

function test_string_serialize_function(M) {
  // Round-trip a Function via StringSerializer (the type-typed path,
  // not just .serialize()).
  const x = M.SX.sym("x");
  const f = M.Function("squared", [x], [M.times(x, x)]);
  const ss = new M.StringSerializer();
  ss.pack(f);
  const blob = ss.encode();
  assertTrue(blob.length > 100, "Function pack blob non-trivial");
  const sd = new M.StringDeserializer(blob);
  const back = sd.unpack_function();
  assertEqual(back.name(), "squared", "unpacked Function name");
  checkArray(back(5), [25], 10, "unpacked Function(5) == 25");
}

function test_string_serialize_sx(M) {
  // Round-trip an SX expression.
  const x = M.SX.sym("x", 3);
  const expr = M.times(x, x);
  const ss = new M.StringSerializer();
  ss.pack(expr);
  const blob = ss.encode();
  const sd = new M.StringDeserializer(blob);
  const back = sd.unpack_sx();
  assertEqual(back.size1(), 3n, "unpacked SX rows");
  assertEqual(back.size2(), 1n, "unpacked SX cols");
}

function test_string_serialize_sparsity(M) {
  // Round-trip a Sparsity.
  const sp = M.Sparsity.lower(5n);
  const ss = new M.StringSerializer();
  ss.pack(sp);
  const blob = ss.encode();
  const sd = new M.StringDeserializer(blob);
  const back = sd.unpack_sparsity();
  assertEqual(back.size1(), 5n, "unpacked Sparsity size1");
  assertEqual(back.size2(), 5n, "unpacked Sparsity size2");
  assertEqual(back.nnz(),  15n, "unpacked Sparsity nnz");
}

function test_unpack_generic_dispatch(M) {
  // StringDeserializer.unpack() should dispatch on pop_type() and route to
  // the right blind_unpack_<type>() method without needing the caller to
  // know the type up-front.  Mirrors the Python `de.unpack()` ergonomics.
  const x = M.SX.sym("x");
  const f = M.Function("squared", [x], [M.times(x, x)]);
  const ss = new M.StringSerializer();
  ss.pack(M.DM([7, 8]));   // -> tag 2  (DM)
  ss.pack(f);              // -> tag 5  (Function)
  ss.pack(M.SX.sym("y", 2)); // -> tag 3  (SX_v1)
  ss.pack(M.Sparsity.dense(2n, 3n)); // -> tag 0 (Sparsity)
  const blob = ss.encode();

  const sd = new M.StringDeserializer(blob);
  const back_dm = sd.unpack();
  assertTrue(back_dm && back_dm.nonzeros, "first unpack() returned DM");
  assertArrayAlmostEqual(back_dm.nonzeros(), [7, 8], 10, "DM via generic unpack");

  const back_fn = sd.unpack();
  assertEqual(back_fn.name(), "squared", "second unpack() returned Function");
  checkArray(back_fn(4), [16], 10, "fn(4) == 16 via generic unpack");

  const back_sx = sd.unpack();
  assertEqual(back_sx.size1(), 2n, "third unpack() returned SX with rows=2");

  const back_sp = sd.unpack();
  assertEqual(back_sp.size1(), 2n, "fourth unpack() returned Sparsity size1");
  assertEqual(back_sp.size2(), 3n, "fourth unpack() returned Sparsity size2");
}

function test_sparsity_serialize_roundtrip(M) {
  // Sparsity.serialize() / .deserialize() round-trip.
  const sp = M.Sparsity.dense(4n, 5n);
  const blob = sp.serialize();
  assertTrue(typeof blob === "string", "Sparsity.serialize -> string");
  const back = M.Sparsity.deserialize(blob);
  assertEqual(back.size1(), 4n, "deserialized size1");
  assertEqual(back.size2(), 5n, "deserialized size2");
  assertEqual(back.nnz(), 20n, "deserialized nnz");
}

runTests([
  ["test_function_serialize_roundtrip", test_function_serialize_roundtrip],
  ["test_string_serialize_dm",          test_string_serialize_dm],
  ["test_string_serialize_function",    test_string_serialize_function],
  ["test_string_serialize_sx",          test_string_serialize_sx],
  ["test_string_serialize_sparsity",    test_string_serialize_sparsity],
  ["test_sparsity_serialize_roundtrip", test_sparsity_serialize_roundtrip],
  ["test_unpack_generic_dispatch",      test_unpack_generic_dispatch],
]);
