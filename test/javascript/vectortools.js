// vectortools.js -- port of test/python/vectortools.py (1 test).
const { assertThrows, checkArray, runTests } = require("./_helpers");

function test_complement(M) {
  // complement([2,1,4,6], 8) -> [0,3,5,7].  Indices out of range error.
  assertThrows(() => M.complement([2, 1, 4, 6], 3n), "complement([...], 3) bounds err");
  assertThrows(() => M.complement([2, 1, 4, 6], 6n), "complement([...], 6) bounds err");
  const wc = M.complement([2, 1, 4, 6], 8n);
  // wc is a vector<casadi_int>; iterate via Array.from / .at(i)
  const arr = Array.from(wc).map(Number);
  // Order may not be guaranteed but should contain the right elements.
  const expected = [0, 3, 5, 7];
  if (arr.length !== expected.length) {
    throw new Error(`complement length: expected ${expected.length}, got ${arr.length}`);
  }
  // Compare as sorted arrays in case order differs.
  const sortedActual = [...arr].sort((a, b) => a - b);
  for (let i = 0; i < expected.length; ++i) {
    if (sortedActual[i] !== expected[i]) {
      throw new Error(`complement[${i}]: expected ${expected[i]}, got ${sortedActual[i]}`);
    }
  }
}

runTests([
  ["test_complement", test_complement],
]);
