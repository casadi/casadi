const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const {spawnSync} = require('node:child_process');
const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'casadi-manifest-'));
try {
  const rules = path.resolve(__dirname, '../../swig/wasm-js/plugins.cmake');
  function configure(extra) {
    const script = path.join(dir, 'test.cmake');
    fs.writeFileSync(script, `cmake_minimum_required(VERSION 3.16)\ninclude("${rules}")\n` +
      'set_property(GLOBAL PROPERTY CASADI_PLUGINS "Nlpsol::fatrop;Importer::shell")\n' +
      extra + `\nwasmjs_write_manifest("${dir}/plugins.json")\n`);
    return spawnSync('cmake', ['-P', script], {encoding: 'utf8'});
  }
  let result = configure('wasmjs_exclude_plugin(importer shell "Needs native compiler")');
  assert.notEqual(result.status, 0);
  assert.match(result.stderr, /nlpsol::fatrop has no Wasm/);
  result = configure('wasmjs_register_plugin(nlpsol fatrop)\nwasmjs_exclude_plugin(importer shell "Needs native compiler")');
  assert.equal(result.status, 0, result.stderr);
  const manifest = JSON.parse(fs.readFileSync(path.join(dir, 'plugins.json')));
  assert.deepEqual(manifest.plugins[0], {kind: 'nlpsol', name: 'fatrop', file: 'libcasadi_nlpsol_fatrop.so'});
  result = configure('wasmjs_register_plugin(nlpsol fatrop)\nwasmjs_exclude_plugin(importer shell [=[Needs "native"; compiler\non host]=])');
  assert.equal(result.status, 0, result.stderr);
  assert.equal(JSON.parse(fs.readFileSync(path.join(dir, 'plugins.json'))).plugins[1].excluded,
    'Needs "native"; compiler\non host');
  result = configure('wasmjs_register_plugin(nlpsol fatrop)\nwasmjs_exclude_plugin(nlpsol fatrop "Conflict")');
  assert.notEqual(result.status, 0);
  assert.match(result.stderr, /both shipped and excluded/);
  console.log('ok -- configure rejects omissions and conflicting exclusions');
} finally {fs.rmSync(dir, {recursive: true, force: true});}
