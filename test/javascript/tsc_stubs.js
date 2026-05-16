// tsc_stubs.js -- TypeScript analog of test/python/pyright_stubs.py.
//
// Runs the system `tsc` compiler over the SWIG-generated casadi.d.ts
// to catch stub-quality regressions: duplicate identifiers, broken
// interface merges, recursive-alias loops, missing `extends` chains,
// wrong overload return-type narrowing under input-coercion unions.
//
// Two tiers, mirroring pyright_stubs.py:
//
//   1. test_dts_compiles        -- tsc against the bare casadi.d.ts.
//                                  Catches the structural defects:
//                                  TS2300 (duplicate identifier),
//                                  TS2339 (missing property),
//                                  TS2456 (circular alias).
//   2. test_types_suite         -- tsc over test/javascript/types/*.ts
//                                  (consumer fixtures that import from
//                                  casadi).  Catches return-type
//                                  narrowing regressions and overload-
//                                  resolution drift.
//
// Both tiers assert tsc emits at most TSC_ERROR_BUDGET errors.
// Ratchet DOWN as stub quality improves; never raise it silently.
//
// Skips gracefully (exit 2) if tsc isn't on PATH, mirroring the
// _have_pyright() guard in pyright_stubs.py.  Skips if casadi.d.ts
// hasn't been built.

const fs = require("fs");
const path = require("path");
const os = require("os");
const cp = require("child_process");

const TSC_ERROR_BUDGET = 0;

const wasmDir = path.resolve(__dirname, "../../build-wasm/swig/wasm-js");
const dtsPath = path.join(wasmDir, "casadi.d.ts");

function haveTsc() {
  const r = cp.spawnSync("tsc", ["--version"], { encoding: "utf8" });
  return r.status === 0;
}

function tscArgs(...extra) {
  return [
    "--noEmit",
    "--strict",
    "--skipLibCheck",
    "--target", "ES2020",
    "--module", "commonjs",
    "--moduleResolution", "node",
    ...extra,
  ];
}

function countErrors(stdout) {
  // tsc lines look like:  file.ts(L,C): error TSxxxx: message
  const matches = stdout.match(/^.+\(\d+,\d+\): error TS\d+:/gm);
  return matches ? matches.length : 0;
}

function runTsc(args, cwd) {
  const r = cp.spawnSync("tsc", args, {
    cwd,
    encoding: "utf8",
    timeout: 120000,
  });
  // tsc emits diagnostics on stdout, not stderr.
  return { status: r.status, stdout: r.stdout || "", stderr: r.stderr || "" };
}

// --- Tier 1: bare casadi.d.ts typechecks under --strict --skipLibCheck ---
function test_dts_compiles() {
  const r = runTsc(tscArgs(dtsPath), os.tmpdir());
  const errors = countErrors(r.stdout);
  if (errors > TSC_ERROR_BUDGET) {
    const sample = r.stdout.split("\n").slice(0, 40).join("\n");
    throw new Error(
      `tsc reported ${errors} error(s) on casadi.d.ts ` +
      `(budget ${TSC_ERROR_BUDGET}). First diagnostics:\n${sample}`);
  }
  if (errors < TSC_ERROR_BUDGET) {
    console.warn(
      `WARN: tsc clean on casadi.d.ts (${errors} errors) -- ` +
      `consider lowering TSC_ERROR_BUDGET in tsc_stubs.js`);
  }
  console.log(`ok 1 - test_dts_compiles  (${errors}/${TSC_ERROR_BUDGET} errors)`);
}

// --- Tier 2: suite over test/javascript/types/*.ts ---
//
// Fixtures `import ... from "./casadi"`, so we copy the d.ts next to
// them in a tempdir and run tsc over the whole dir.  This avoids any
// tsconfig/path-mapping dance and keeps fixtures resolution-agnostic.
function test_types_suite() {
  const typesDir = path.join(__dirname, "types");
  if (!fs.existsSync(typesDir)) {
    console.log(`ok 2 - test_types_suite  (skip: no types/ dir)`);
    return;
  }
  const fixtures = fs.readdirSync(typesDir).filter((f) => f.endsWith(".ts"));
  if (fixtures.length === 0) {
    console.log(`ok 2 - test_types_suite  (skip: no .ts fixtures)`);
    return;
  }
  const tmp = fs.mkdtempSync(path.join(os.tmpdir(), "casadi-tsc-"));
  try {
    fs.copyFileSync(dtsPath, path.join(tmp, "casadi.d.ts"));
    for (const f of fixtures) {
      fs.copyFileSync(path.join(typesDir, f), path.join(tmp, f));
    }
    const r = runTsc(tscArgs("--noUnusedLocals", "false", ...fixtures.map((f) => f)), tmp);
    const errors = countErrors(r.stdout);
    if (errors > TSC_ERROR_BUDGET) {
      // Rewrite the tempdir prefix back to test/javascript/types/ in
      // the diagnostic so failures point at the source fixtures.
      const remapped = r.stdout.replaceAll(tmp + path.sep, "test/javascript/types/");
      const sample = remapped.split("\n").slice(0, 40).join("\n");
      throw new Error(
        `tsc reported ${errors} error(s) over types/ fixtures ` +
        `(budget ${TSC_ERROR_BUDGET}). First diagnostics:\n${sample}`);
    }
    if (errors < TSC_ERROR_BUDGET) {
      console.warn(
        `WARN: tsc clean over types/ (${errors} errors) -- ` +
        `consider lowering TSC_ERROR_BUDGET in tsc_stubs.js`);
    }
    console.log(`ok 2 - test_types_suite  (${fixtures.length} fixtures, ${errors}/${TSC_ERROR_BUDGET} errors)`);
  } finally {
    fs.rmSync(tmp, { recursive: true, force: true });
  }
}

// --- Main ---
(function main() {
  if (!haveTsc()) {
    console.error("SKIP: tsc not on PATH (install with `npm i -g typescript`)");
    process.exit(2);
  }
  if (!fs.existsSync(dtsPath)) {
    console.error(`SKIP: ${dtsPath} not built. ` +
      `Run cmake --build build-wasm --target casadi_wasm first.`);
    process.exit(2);
  }
  try {
    test_dts_compiles();
    test_types_suite();
    console.log("1..2");
  } catch (e) {
    console.error(`not ok - ${e.message}`);
    process.exit(1);
  }
})();
