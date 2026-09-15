  // Fetch each plugin once, including concurrent loads; permit retries after failures.
  const __pluginLoads = new Map();
  const __loadPlugin = (kind, name, register) => {
    if (!/^[A-Za-z0-9_]+$/.test(name)) {
      return Promise.reject(new Error("Invalid plugin name: " + name));
    }
    const soname = "libcasadi_" + kind + "_" + name + ".so";
    if (!__pluginLoads.has(soname)) {
      const pending = Promise.resolve().then(async () => {
        if (!M.FS.analyzePath("/" + soname).exists) {
          let bytes;
          if (typeof process !== "undefined" && process.versions && process.versions.node) {
            bytes = require("fs").readFileSync(__path.join(__dirname, soname));
          } else {
            const url = M.locateFile ? M.locateFile(soname, "") : soname;
            const response = await fetch(url);
            if (!response.ok) throw new Error("Failed to fetch plugin " + soname + ": " + response.status);
            bytes = new Uint8Array(await response.arrayBuffer());
          }
          M.FS.writeFile("/" + soname, bytes);
        }
        return register(name);
      });
      __pluginLoads.set(soname, pending);
      pending.catch(() => {
        __pluginLoads.delete(soname);
        // A failed registration may have left an unusable file behind.
        try { M.FS.unlink("/" + soname); } catch (_) {}
      });
    }
    return __pluginLoads.get(soname);
  };

  for (const kind of ["nlpsol", "conic", "linsol", "integrator", "rootfinder",
      "interpolant", "expm", "dple", "blas", "filesystem"]) {
    const fn = "load_" + kind;
    const original = __m[fn];
    if (typeof original === "function") {
      __m[fn] = name => __loadPlugin(kind, name, original.bind(__m));
    }
  }
  // Class-based entry points need the same asynchronous fetch as free functions.
  for (const [cls, kind, free] of [
      ["Linsol", "linsol", "load_linsol"], ["XmlFile", "xmlfile"],
      ["Importer", "importer"], ["ModelicaParser", "modelicaparser"]]) {
    if (!__m[cls] || typeof __m[cls].load_plugin !== "function") continue;
    const original = __m[cls].load_plugin.bind(__m[cls]);
    __m[cls].load_plugin = free && __m[free]
      ? __m[free] : name => __loadPlugin(kind, name, original);
  }
