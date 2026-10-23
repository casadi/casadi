#!/usr/bin/env python
"""Pair non-Python \\doctest blocks of the C++ headers with the Python block they translate.

Writes JSON {normalized python block: [[lang, block text], ...]} for the Sphinx hook, which inserts
those blocks after the matching Python example on the API page. A block in another language belongs
to the closest preceding Python block of the same Doxygen comment.
Usage: python doctest_sidecar.py <casadi source dir> <output.json>
"""
import json
import os
import re
import sys


def norm(text):
    return "\n".join(l.strip() for l in text.strip().split("\n") if l.strip())


def pairs(path):
    out = {}
    lines = open(path).read().split("\n")
    in_comment, cur, last_py = False, None, None
    for line in lines:
        if "/**" in line:
            in_comment, last_py = True, None
        s = re.sub(r"^\s*\*?\s?", "", line)
        m = re.match(r"\\doctest(?:\{(\w+)\})?\s*$", s.strip())
        if in_comment and m:
            cur = (m.group(1) or "python", [])
        elif cur is not None and re.sub(r"\*/\s*$", "", s).strip() == "\\enddoctest":
            lang, body = cur[0], "\n".join(cur[1]).strip("\n")
            if lang == "python":
                last_py = norm(body)
            elif last_py is not None:
                out.setdefault(last_py, []).append([lang, body])
            cur = None
        elif cur is not None:
            cur[1].append(s)
        if "*/" in line:
            in_comment, cur = False, None
    return out


if __name__ == "__main__":
    root, dest = sys.argv[1], sys.argv[2]
    data = {}
    for dirpath, _, files in os.walk(root):
        for f in sorted(files):
            if f.endswith(".hpp"):
                for k, v in pairs(os.path.join(dirpath, f)).items():
                    data.setdefault(k, []).extend(v)
    json.dump(data, open(dest, "w"), indent=0)
    print("doctest_sidecar: %d python blocks with translations" % len(data))
