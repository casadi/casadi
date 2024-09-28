import os
import sys
import subprocess

identity = sys.argv[1]
dir = sys.argv[2]

# Recursively look for all shared libraries in `dir`
for root, dirs, files in os.walk(dir):
    for file in files:
        if file.endswith(".dylib") or file.endswith(".mexmaci64") or file.endswith(".mexmaca64") or file.endswith(".so") or file.endswith(".mex") or os.access(file, os.X_OK):
            # Get the full path to the library
            path = os.path.join(root, file)

            # codesign the library
            subprocess.call(["codesign", "--remove-signature", path])
            subprocess.call(["codesign", "--force", "--sign", identity, path])
