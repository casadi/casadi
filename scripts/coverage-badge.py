#!/usr/bin/env python3
"""
Script that extracts the test coverage percentage and saves it as a JSON
shields.io endpoint.
"""

from os.path import join, normpath, dirname, realpath
import re

script_dir = dirname(realpath(__file__))
cov_dir = normpath(join(dirname(script_dir), "docs", "Coverage"))

json = """\
{{
  "schemaVersion": 1,
  "label": "Test Coverage",
  "message": "{linecov}%",
  "color": "green"
}}
"""

def main():
    with open(join(cov_dir, 'index.html'), 'r') as f:
        pattern = r'<td class="headerCovTableEntry\w+">([\d.]+)'
        linecov, = map(lambda m: m.group(1),
                       re.finditer(pattern, f.read()))
        print(linecov)

    with open(join(cov_dir, "shield.io.coverage.json"), 'w') as f:
        f.write(json.format(linecov=linecov))


if __name__ == "__main__":
    main()
