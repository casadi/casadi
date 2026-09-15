"""Stage exactly the plugins declared by the Wasm build, excluding stale files."""
import json
from pathlib import Path
import shutil
import sys


def stage(source, destination):
    if source.resolve() == destination.resolve():
        raise ValueError("Source and destination must be different directories")
    manifest = source / "plugins.json"
    plugins = json.loads(manifest.read_text())["plugins"]
    files = [source / p["file"] for p in plugins if "file" in p]
    for file in files:
        if not file.is_file():
            raise FileNotFoundError(f"Required Wasm plugin is missing: {file}")
    destination.mkdir(parents=True, exist_ok=True)
    for file in destination.glob("libcasadi_*.so"):
        file.unlink()
    for file in files:
        shutil.copy2(file, destination / file.name)
    shutil.copy2(manifest, destination / manifest.name)


if __name__ == "__main__":
    stage(Path(sys.argv[1]), Path(sys.argv[2]))
