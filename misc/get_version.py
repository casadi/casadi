version_parts = {}

import sys
import os
import re

cd = os.getcwd()

if len(sys.argv)>1:
  cd = sys.argv[1]

with open(os.path.join(cd,"CMakeLists.txt"),"r") as f:
  for line in f.readlines():
    if "_VERSION" in line:
      m = re.search(r"set\(CASADI_(\w+)_VERSION (\d+)", line)
      if m:
        version_parts[m.group(1).lower()] = m.group(2)

sys.stdout.write("{major}.{minor}.{patch}".format(**version_parts))
