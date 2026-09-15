// Run with Playwright installed; use two origins to exercise CDN-style loading.
const {chromium} = require('playwright');
const http = require('node:http');
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');
const dir = path.resolve(process.argv[2] || path.join(__dirname, '../../build-wasm/swig/wasm-js'));
const requests = [];
function serve(handler) {
  const server = http.createServer(handler);
  return new Promise(resolve => server.listen(0, '127.0.0.1', () => resolve(server)));
}
(async () => {
  let browser, assets, pages, timer;
  try {
    assets = await serve((req, res) => {
      requests.push(req.url);
      const name = req.url.replace(/^\/package\//, '');
      if (!req.url.startsWith('/package/') || name !== path.basename(name)) {
        res.writeHead(404); res.end(); return;
      }
      const file = path.join(dir, name);
      if (!fs.existsSync(file)) {res.writeHead(404); res.end(); return;}
      res.setHeader('Access-Control-Allow-Origin', '*');
      res.setHeader('Content-Type', name.endsWith('.wasm') || name.endsWith('.so')
        ? 'application/wasm' : name.endsWith('.json') ? 'application/json' : 'text/javascript');
      fs.createReadStream(file).pipe(res);
    });
    pages = await serve((req, res) => {
      if (req.url === '/loader.js' || req.url === '/checks.js') {
        res.setHeader('Content-Type', 'text/javascript');
        fs.createReadStream(req.url === '/loader.js'
          ? path.join(__dirname, '../../docs/examples/javascript/_casadi_browser.js')
          : path.join(__dirname, 'wasm_plugins.js')).pipe(res);
      } else {
        res.setHeader('Content-Type', 'text/html');
        res.end('<script src="/loader.js"></script><script src="/checks.js"></script>');
      }
    });
    browser = await chromium.launch({headless: true});
    const page = await browser.newPage();
    const errors = [];
    page.on('pageerror', e => errors.push(String(e)));
    await page.goto(`http://127.0.0.1:${pages.address().port}/nested/page.html`);
    const base = `http://127.0.0.1:${assets.address().port}/package/`;
    const result = await Promise.race([page.evaluate(async base => {
      const ca = await loadCasadi(base);
      const manifest = await fetch(base + 'plugins.json').then(r => r.json());
      const loaded = await testWasmPlugins(ca, manifest);
      await Promise.all([ca.Linsol.load_plugin('qr'), ca.load_linsol('qr')]);
      return loaded;
    }, base), new Promise((_, reject) => {
      timer = setTimeout(() => reject(Error('Browser plugin checks timed out')), 120000);
    })]);
    assert.deepEqual(errors, []);
    assert(result.length > 0);
    assert(requests.includes('/package/libcasadi_linsol_qr.so'));
    assert.equal(requests.filter(p => p === '/package/libcasadi_linsol_qr.so').length, 1);
    console.log('Browser plugin checks passed:', result);
  } finally {
    clearTimeout(timer);
    if (browser) await browser.close();
    if (pages) await new Promise(resolve => pages.close(resolve));
    if (assets) await new Promise(resolve => assets.close(resolve));
  }
})().catch(e => {console.error(e); process.exitCode = 1;});
