import os
import sys
import subprocess

identity = sys.argv[1]
dir = sys.argv[2]

# Recursively look for all shared libraries in `dir`
for root, dirs, files in os.walk(dir):
    for file in files:
        if file.endswith(".dylib"):
            # Get the full path to the library
            path = os.path.join(root, file)

            # codesign the library
            subprocess.call(["codesign", "--force", "--sign", identity, path])
