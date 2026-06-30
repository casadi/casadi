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

// --- Tier 3: typecheck the .js runtime suite via --allowJs --checkJs ---
//
// Each .js test loads the casadi module via `require(modulePath)()`
// where modulePath is computed dynamically.  tsc can't statically
// resolve dynamic require args, so we wrap each .js in a tempdir
// copy with a `// @ts-check` directive prepended AND a JSDoc-typed
// shim that aliases `M` to the correct module shape.  Originals
// stay untouched.
//
// The wrapper rewrites every `const M = await create();` line to
//   /** @type {Awaited<ReturnType<typeof import("./casadi").default>>} */
//   const M = await create();
// so M-rooted access (`M.DM.sym(...)` etc.) gets full type info.
//
// Files explicitly opted-out (in JS_TYPECHECK_SKIP) are left alone.
// Helper / non-test files (filename starting with `_`) are also
// skipped: they have no `M = await create()` line to annotate.
//
// Per-target ratchet via JS_TIER3_BUDGET because the .js suite was
// not written with TS-correctness in mind -- expect a non-zero
// initial count.  Ratchet DOWN as files get cleaned up.
// Tier-3 budget is non-zero because expanding to runTests-style files
// (integration.js, sx.js, matrix.js, etc.) surfaced ~44 d.ts↔runtime
// mismatches that need separate cleanup:
//
//   * `.T` is a getter at runtime (casadi.js %insert("js") block) but
//     the d.ts declares `T(): DM` as a method.  Tests do `a.T` not `a.T()`.
//   * StringDeserializer's `unpack` is added dynamically at runtime
//     (JS dispatcher on pop_type()).  Not in the d.ts.
//   * Custom GenericType field access (e.g. `stats.bar`) needs narrowing
//     in the test code; could be improved with type guards.
//   * A couple of test-code bugs: passing function references instead of
//     calling them, treating returned arrays as scalars.
//   * Director subclass internal `_ptr` access.
//
// All real issues, none structural to the d.ts shape.  Ratchet down as
// each category is addressed; never raise silently.
const JS_TIER3_BUDGET = 44;
const JS_TYPECHECK_SKIP = new Set([
  "tsc_stubs.js",      // this harness itself
  "_helpers.js",       // helper module
  "_leak_check.js",    // requires --expose-gc
]);

function getNodeTypeRoots() {
  // Locate @types/node on the system.  Debian/Ubuntu uses
  // /usr/share/nodejs/@types; node_modules/@types is the usual local
  // path.  Fall back to skipping if neither is found.
  const candidates = [
    "/usr/share/nodejs/@types",
    path.join(process.cwd(), "node_modules/@types"),
  ];
  for (const c of candidates) {
    if (fs.existsSync(path.join(c, "node"))) return c;
  }
  return null;
}

function test_js_suite_typecheck() {
  const jsDir = __dirname;
  const all = fs.readdirSync(jsDir).filter((f) =>
    f.endsWith(".js") && !JS_TYPECHECK_SKIP.has(f));
  if (all.length === 0) {
    console.log(`ok 3 - test_js_suite_typecheck  (skip: no .js tests)`);
    return;
  }
  const typeRoots = getNodeTypeRoots();
  if (!typeRoots) {
    console.log(`ok 3 - test_js_suite_typecheck  (skip: no @types/node found)`);
    return;
  }
  const tmp = fs.mkdtempSync(path.join(os.tmpdir(), "casadi-jstsc-"));
  let n_checked = 0;
  const toCheck = [];
  try {
    fs.copyFileSync(dtsPath, path.join(tmp, "casadi.d.ts"));
    // Stub casadi.js so the `require("./casadi")` in our annotated
    // wrapper resolves -- TS only needs the .d.ts to type things;
    // the actual runtime never gets executed.
    fs.writeFileSync(path.join(tmp, "casadi.js"),
      "module.exports = function() { return Promise.resolve({}); };\n");
    // Copy ALL .js files (including helpers / harness skip-list ones)
    // so local relative-requires like `require("./_helpers")` resolve.
    // Only the non-skipped files are passed to tsc for checking.
    for (const f of fs.readdirSync(jsDir).filter((g) => g.endsWith(".js"))) {
      fs.copyFileSync(path.join(jsDir, f), path.join(tmp, f));
    }
    for (const f of all) {
      const src = fs.readFileSync(path.join(jsDir, f), "utf8");
      // Replace the dynamic `require(modulePath)` with a static
      // `require("./casadi")` so TS picks up our local d.ts; and
      // annotate `M = await create()` with the typed-shape JSDoc.
      let body = src;
      // Annotate `create` as the factory function regardless of what
      // the dynamic require() actually resolves to.  The d.ts shape
      // `typeof import("./casadi")` is the module namespace; the
      // factory returns a Promise<that-namespace>.
      body = body.replace(/(const\s+create\s*=\s*require\(modulePath\)\s*;)/,
        `/** @type {() => Promise<typeof import("./casadi")>} */\n  // @ts-ignore -- runtime require is the factory; d.ts has classes as named exports\n  $1`);
      // The .js suite uses a recurring test-runner pattern:
      //   const tests = [["test_a", testA], ...];
      //   for (const [name, fn] of tests) { await fn(M, ...); }
      // Without an annotation, TS infers tests as
      // `(string | testfn)[][]` and the destructured `fn` becomes
      // `string | testfn`, which isn't callable.  This isn't a real
      // type-error in the bindings; it's a JS test-harness shape gap.
      // Annotate `tests` to preserve the tuple shape so iteration
      // typechecks.
      body = body.replace(/(const\s+tests\s*=\s*\[)/,
        `/** @type {Array<[string, (...args: any[]) => any]>} */\n  $1`);
      // Same fix for the plugin / solver dispatcher tables used in
      // conic.js / linsol.js / rootfinder.js -- they're 3-tuples
      // `[plugin_name, opts_factory_or_obj, caps]` that TS-widens to
      // a useless union without annotation.
      body = body.replace(/(const\s+(ALL_PLUGINS|ALL_SOLVERS|ALL_FUNCTIONS)\s*=\s*\[)/,
        `/** @type {any[]} */\n$1`);
      // Files that delegate test-running to _helpers.runTests have no
      // local `const M = await create()` -- M arrives via the test
      // function's first parameter.  Annotate `function test_XXX(M, ...)`
      // declarations so M is typed against the module shape.  The
      // d.ts now emits classes in callable+newable form (see
      // wasm_js.cxx `class <Name>__class` + `const <Name>` pattern),
      // so the JS-test idiom `M.DM(x)` (sans `new`) typechecks.
      body = body.replace(/^function\s+(test_\w+)\s*\(\s*M([,)])/gm,
        `/** @param {typeof import("./casadi")} M */\nfunction $1(M$2`);
      // Skip files that don't have the create() pattern (likely
      // helper modules that already passed the _ prefix filter).
      // Skip files where no annotation pattern matched -- nothing to
      // typecheck meaningfully.  Includes helper / smoke files that
      // load casadi via inline paths or don't take M as a parameter.
      if (body === src) continue;
      n_checked++;
      // Overwrite the plain copy with the annotated version.
      fs.writeFileSync(path.join(tmp, f), `// @ts-check\n${body}`);
      toCheck.push(f);
    }
    if (n_checked === 0) {
      console.log(`ok 3 - test_js_suite_typecheck  (skip: no annotatable .js files)`);
      return;
    }
    // Pass ONLY the annotated files to tsc (so helpers stay typecheck-
    // exempt but resolve via local require paths).
    const args = [
      "--noEmit", "--allowJs", "--checkJs", "--skipLibCheck",
      "--target", "ES2020", "--module", "commonjs",
      "--moduleResolution", "node",
      "--typeRoots", typeRoots, "--types", "node",
      ...toCheck,
    ];
    const r = runTsc(args, tmp);
    const errors = countErrors(r.stdout);
    if (errors > JS_TIER3_BUDGET) {
      const remapped = r.stdout.replaceAll(tmp + path.sep, "test/javascript/");
      const sample = remapped.split("\n").slice(0, 60).join("\n");
      throw new Error(
        `tsc --checkJs reported ${errors} error(s) over .js suite ` +
        `(${n_checked} files, budget ${JS_TIER3_BUDGET}). First diagnostics:\n${sample}`);
    }
    if (errors < JS_TIER3_BUDGET) {
      console.warn(
        `WARN: tsc --checkJs clean (${errors} errors) -- ` +
        `consider lowering JS_TIER3_BUDGET`);
    }
    console.log(`ok 3 - test_js_suite_typecheck  (${n_checked} files, ${errors}/${JS_TIER3_BUDGET} errors)`);
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
    test_js_suite_typecheck();
    console.log("1..3");
  } catch (e) {
    console.error(`not ok - ${e.message}`);
    process.exit(1);
  }
})();
