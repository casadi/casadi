"""Verify that a transitive solver dependency invalidates the Wasm side module."""
from pathlib import Path
import subprocess
import sys

build = Path(sys.argv[1]).resolve()
side = build / "swig/wasm-js/libcasadi_nlpsol_fatrop.so"
leaf = build / "external_projects/lib/libblasfeo.a"
before = side.stat().st_mtime_ns
if not leaf.is_file():
    raise FileNotFoundError(leaf)
leaf.touch()
result = subprocess.run(
    ["cmake", "--build", str(build), "--target", "casadi_wasm", "-j2"],
    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
if result.returncode:
    print(result.stdout)
    sys.exit(result.returncode)
if side.stat().st_mtime_ns <= before:
    raise AssertionError("Updating BLASFEO did not relink the Fatrop side module")
print("ok -- updating BLASFEO relinks Fatrop through its transitive dependency")
