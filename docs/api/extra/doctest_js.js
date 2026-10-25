// Run JavaScript doctest blocks (node REPL style: `> expr` + expected output) against casadi-wasm.
//
// Usage: node docs/api/extra/doctest_js.js <blocks.txt>
// The file holds blocks separated by lines `#### <name>`; each block is the text between
// \doctest and \enddoctest of a \lang{javascript} group. The module is build-wasm/swig/wasm-js/casadi.js
// (override with CASADI_WASM_JS). Exit code = number of failing blocks.
const fs = require("fs");
const path = require("path");
const util = require("util");
const vm = require("vm");

const modulePath = process.env.CASADI_WASM_JS ||
  path.resolve(__dirname, "..", "..", "..", "build-wasm", "swig", "wasm-js", "casadi.js");

function parseBlocks(text) {
  const blocks = [];
  let name = null, buf = [];
  for (const line of text.split("\n")) {
    if (line.startsWith("#### ")) {
      if (name !== null) blocks.push([name, buf.join("\n")]);
      name = line.slice(5).trim(); buf = [];
    } else {
      buf.push(line);
    }
  }
  if (name !== null) blocks.push([name, buf.join("\n")]);
  return blocks;
}

// `> expr` starts an example, `... ` continues it, following non-blank lines are its output.
function parseExamples(text) {
  const out = [];
  const lines = text.split("\n");
  let i = 0;
  while (i < lines.length) {
    if (lines[i].startsWith("> ")) {
      const src = [lines[i].slice(2)]; i++;
      while (i < lines.length && lines[i].startsWith("... ")) { src.push(lines[i].slice(4)); i++; }
      const exp = [];
      while (i < lines.length && !lines[i].startsWith("> ") && lines[i].trim() !== "") { exp.push(lines[i]); i++; }
      out.push([src.join("\n"), exp.join("\n")]);
    } else {
      i++;
    }
  }
  return out;
}

const norm = s => s.split("\n").map(l => l.trim()).filter(l => l).join("\n");

// Numbers as in doctest_headers.py: half a unit in the last written digit, at least 1e-9 relative
const NUMBER = /(?<![\w.@])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?(?![\w.])/g;
const numbers = s => [...s.matchAll(NUMBER)].filter(m => !/^[-+]?0\d/.test(m[0]));

function tolerance(lit) {
  const m = lit.match(/^[-+]?\d*(?:\.(\d*))?(?:[eE]([-+]?\d+))?$/);
  const tol = m[1] === undefined && m[2] === undefined ? 0 : 0.5 * 10 ** (Number(m[2] || 0) - (m[1] || "").length);
  return Math.max(tol, 1e-9 * (Math.abs(Number(lit)) || 1));
}

function asExpected(expected, got) {
  const w = numbers(expected), g = numbers(got);
  if (w.length !== g.length) return got;
  let out = "", pos = 0;
  w.forEach((a, k) => {
    const b = g[k];
    if (Math.abs(Number(a[0]) - Number(b[0])) <= tolerance(a[0]) * (1 + 1e-9)) {
      out += got.slice(pos, b.index) + a[0];
      pos = b.index + b[0].length;
    }
  });
  return out + got.slice(pos);
}

function matches(expected, got) {
  const e = norm(expected), g = asExpected(e, norm(got));
  if (e === g) return true;
  if (!e.includes("...")) return false;
  const pat = e.split("...").map(x => x.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")).join("[\\s\\S]*");
  return new RegExp("^" + pat + "$").test(g);
}

async function main() {
  const create = require(modulePath);
  const ca = await create();
  const blocks = parseBlocks(fs.readFileSync(process.argv[2], "utf8"));
  let nfail = 0;
  for (const [name, text] of blocks) {
    let printed = "";
    const context = vm.createContext({
      ca, casadi: ca,
      console: { log: (...a) => { printed += a.map(x => typeof x === "string" ? x : util.inspect(x)).join(" ") + "\n"; } },
    });
    let ok = true;
    for (const [src, expected] of parseExamples(text)) {
      printed = "";
      let got;
      try {
        const result = vm.runInContext(src, context, { filename: name });
        got = printed + (result === undefined ? "" : util.inspect(result));
      } catch (err) {
        got = printed + "Uncaught " + String(err);
      }
      if (!matches(expected, got)) {
        ok = false;
        console.log("**********************************************************************");
        console.log(`${name}\nFailed example:\n    ${src.replace(/\n/g, "\n    ")}\nExpected:\n    ${expected.replace(/\n/g, "\n    ")}\nGot:\n    ${got.trimEnd().replace(/\n/g, "\n    ")}`);
      }
    }
    if (!ok) nfail++;
  }
  console.log(`doctest_js: ${blocks.length} blocks, ${nfail} failures`);
  process.exit(nfail);
}

main();
