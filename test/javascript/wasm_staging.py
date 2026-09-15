"""Regression checks for the release staging contract."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location(
    "stage_plugins", Path(__file__).resolve().parents[2] / "swig/wasm-js/stage_plugins.py")
staging = importlib.util.module_from_spec(spec)
spec.loader.exec_module(staging)


class StagingTests(unittest.TestCase):
    def test_manifest_controls_distribution(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, destination = Path(tmp) / "build", Path(tmp) / "dist"
            source.mkdir()
            destination.mkdir()
            name = "libcasadi_nlpsol_fatrop.so"
            (source / name).write_bytes(b"wasm")
            (source / "libcasadi_conic_old.so").write_bytes(b"stale build")
            (destination / "libcasadi_conic_old.so").write_bytes(b"stale dist")
            (source / "plugins.json").write_text(json.dumps({"plugins": [
                {"kind": "nlpsol", "name": "fatrop", "file": name},
                {"kind": "importer", "name": "shell", "excluded": "Native compiler"}]}))
            staging.stage(source, destination)
            self.assertEqual(sorted(p.name for p in destination.iterdir()),
                             [name, "plugins.json"])
            (source / name).unlink()
            with self.assertRaises(FileNotFoundError):
                staging.stage(source, destination)
            self.assertEqual((destination / name).read_bytes(), b"wasm")


if __name__ == "__main__":
    unittest.main()
