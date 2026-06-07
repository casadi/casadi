//
//     MIT No Attribution
//
//     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//     Permission is hereby granted, free of charge, to any person obtaining a copy of this
//     software and associated documentation files (the "Software"), to deal in the Software
//     without restriction, including without limitation the rights to use, copy, modify,
//     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//     permit persons to whom the Software is furnished to do so.
//
//     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
// Browser loader for the SWIG -wasm-js `casadi.js` wrapper.
//
// `casadi.js` is generated for Node (CommonJS): its first lines do
//   const createWasm = require('./casadi_wasm.js');
//   const __path     = require('path');
//   module.exports   = async function createcasadi() { ... __dirname ... }
// ...and the remaining ~25k lines are plain ECMAScript.  Rather than
// patch the generator, we evaluate it here inside a tiny CommonJS
// sandbox that supplies exactly those three node-isms, so the same
// unmodified file runs in a browser.
//
// Why this is NOT a CORS problem: everything is fetched same-origin
// from the page's own server.  The only requirement is that the page
// is *served over http(s)*, not opened as a file:// URL -- browsers
// refuse to fetch() a sibling .wasm over file://.  See README.md.
//
// Usage (see _template.html):
//   <script src="_casadi_browser.js"></script>
//   <script src="rosenbrock.js"></script>
//   <script>loadCasadi("./").then(M => example(M, console.log));</script>

async function loadCasadi(base) {
  base = base || "./";
  if (!base.endsWith("/")) base += "/";

  // Most common mistake: opening the page as a file:// URL.  Browsers
  // reject fetch() of local files with an opaque "NetworkError", so
  // detect it up front and explain.
  if (typeof location !== "undefined" && location.protocol === "file:") {
    throw new Error(
      "this page must be served over http(s), not opened as a file:// URL " +
      "(browsers block fetch() of local files).  From this folder run: " +
      "`python3 -m http.server 8000` then open http://localhost:8000/");
  }

  // Evaluate a CommonJS module fetched from `base` inside a sandbox.
  const evalCjs = async (file, requireFn) => {
    let resp;
    try {
      resp = await fetch(base + file);
    } catch (e) {
      throw new Error(
        `could not fetch ${base + file} (${e.message || e}).  Is the page ` +
        `served over http, and are the casadi build artifacts (casadi.js, ` +
        `casadi_wasm.js, casadi_wasm.wasm, casadi_wasm.data) present at ` +
        `"${base}"?  See README.md.`);
    }
    if (!resp.ok) throw new Error(
      `${base + file}: ${resp.status} ${resp.statusText} -- artifact missing ` +
      `at "${base}"?  Symlink/copy them in or pass ?casadi=<url>.  See README.md.`);
    const src = await resp.text();
    const factory = new Function(
      "module", "exports", "require", "__dirname", "__filename",
      src + "\n;return module.exports;");
    const module = { exports: {} };
    return factory(module, module.exports, requireFn, base.replace(/\/$/, ""), base + file);
  };

  // 1) The emscripten core.  It auto-detects the browser environment
  //    (no `process`) and will fetch casadi_wasm.{wasm,data} via the
  //    locateFile we pass from the wrapper -- so it never touches the
  //    `require` we hand it; we still provide one that errors loudly.
  const createWasm = await evalCjs("casadi_wasm.js",
    (p) => { throw new Error("casadi_wasm.js required unexpected module: " + p); });

  // 2) The SWIG wrapper.  Resolve its two requires; `path` only needs
  //    `join` (used by the wrapper's locateFile).
  const requireForWrapper = (p) => {
    if (p === "./casadi_wasm.js" || p.endsWith("/casadi_wasm.js")) return createWasm;
    if (p === "path") return { join: (...a) => a.filter(Boolean).join("/").replace(/\/{2,}/g, "/") };
    throw new Error("casadi.js required unexpected module: " + p);
  };
  const createcasadi = await evalCjs("casadi.js", requireForWrapper);

  return await createcasadi();
}

// Also export for bundlers / module contexts that load this file via require/import.
if (typeof module !== "undefined" && module.exports) module.exports = { loadCasadi };
