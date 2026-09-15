// Exercise the actual injected loader in a browser-like context, without a Wasm build.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../../swig/wasm-js/load_plugins.js'), 'utf8');
(async () => {
  const files = new Map();
  const registrations = [], requests = [];
  let fail = true;
  const ca = {XmlFile: {}, Linsol: {}, Importer: {}, ModelicaParser: {}};
  for (const kind of ['nlpsol', 'conic', 'linsol', 'blas', 'filesystem']) {
    ca['load_' + kind] = name => registrations.push(kind + ':' + name);
  }
  for (const cls of ['XmlFile', 'Linsol', 'Importer', 'ModelicaParser']) {
    ca[cls].load_plugin = name => registrations.push(cls + ':' + name);
  }
  vm.runInNewContext(source, {__m: ca, Uint8Array, M: {
    locateFile: name => 'https://cdn.example/package/' + name,
    FS: {analyzePath: p => ({exists: files.has(p)}),
      writeFile: (p, b) => files.set(p, b), unlink: p => files.delete(p)}
  }, fetch: async url => {
    requests.push(url);
    return {ok: !fail, status: 404, arrayBuffer: async () => new ArrayBuffer(8)};
  }});
  await assert.rejects(ca.load_nlpsol('fatrop'), /404/);
  fail = false;
  await Promise.all([ca.load_nlpsol('fatrop'), ca.load_nlpsol('fatrop')]);
  assert.equal(requests.length, 2);
  assert.deepEqual(registrations, ['nlpsol:fatrop']);
  await ca.load_nlpsol('fatrop');
  assert.equal(requests.length, 2);
  await Promise.all([ca.Linsol.load_plugin('qr'), ca.load_linsol('qr')]);
  assert.equal(registrations.filter(x => x === 'linsol:qr').length, 1);
  await ca.XmlFile.load_plugin('tinyxml');
  await ca.load_blas('blasfeo');
  await ca.load_filesystem('ghc');
  assert(requests.every(url => url.startsWith('https://cdn.example/package/libcasadi_')));
  await assert.rejects(ca.load_nlpsol('../fatrop'), /Invalid plugin name/);
  console.log('ok -- async loading, retries, deduplication, class loaders and CDN URLs');
})().catch(e => {console.error(e); process.exitCode = 1;});

// The example loader must preserve URL schemes and resolve relative asset bases.
(async () => {
  const browserSource = fs.readFileSync(path.join(__dirname,
    '../../docs/examples/javascript/_casadi_browser.js'), 'utf8');
  for (const [base, expected] of [
    ['https://cdn.example/pkg/', 'https://cdn.example/pkg/'],
    ['../assets/', 'https://app.example/assets/']]) {
    const urls = [];
    const context = vm.createContext({URL, location: {
      href: 'https://app.example/nested/page.html', protocol: 'https:'
    }, fetch: async url => {
      urls.push(url);
      return {ok: true, text: async () => url.endsWith('casadi_wasm.js')
        ? 'module.exports = async options => options.locateFile("casadi_wasm.wasm", "");'
        : 'const wasm = require("./casadi_wasm.js"); const path = require("path"); module.exports = () => wasm({locateFile: p => path.join(__dirname, p)});'};
    }});
    vm.runInContext(browserSource, context);
    assert.equal(await context.loadCasadi(base), expected + 'casadi_wasm.wasm');
    assert.deepEqual(urls, [expected + 'casadi_wasm.js', expected + 'casadi.js']);
  }
  console.log('ok -- browser adapter preserves HTTPS and relative asset paths');
})().catch(e => {console.error(e); process.exitCode = 1;});
