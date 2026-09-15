const fs = require('fs');
const path = require('path');
const {testWasmPlugins} = require('./wasm_plugins');
(async () => {
  const dir = path.resolve(process.argv[2] || path.join(__dirname, '../../build-wasm/swig/wasm-js'));
  const manifest = JSON.parse(fs.readFileSync(path.join(dir, 'plugins.json'), 'utf8'));
  for (const plugin of manifest.plugins) {
    if (!plugin.excluded && !fs.existsSync(path.join(dir, plugin.file))) {
      throw Error('Missing required plugin artifact: ' + plugin.file);
    }
  }
  const ca = await require(path.join(dir, 'casadi.js'))();
  console.log('Verified plugins:', await testWasmPlugins(ca, manifest));
})().catch(e => { console.error(e); process.exitCode = 1; });
