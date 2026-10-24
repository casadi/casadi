#!/usr/bin/env python
"""Pair the non-Python \\lang groups of the C++ headers with the Python example they translate.

Writes JSON {normalized python block: [[lang, parts], ...]} for the Sphinx hook, which inserts
those groups after the matching Python example on the API page. parts is a list of
["text", paragraph text], ["verbatim", text] and ["doctest", text]. A group in another language
belongs to the closest preceding Python group of the same Doxygen comment and follows its last
\\doctest; groups without such an anchor are left out.
Usage: python doctest_sidecar.py <casadi source dir> <output.json>
"""
import json
import os
import sys

from doctest_headers import parse


def norm(text):
    return "\n".join(l.strip() for l in text.strip().split("\n") if l.strip())


def pairs(path):
    out = {}
    comment, last_py = None, None
    for g in parse(path):
        if g["comment"] != comment:
            comment, last_py = g["comment"], None
        if g["lang"] == "python":
            last_py = None
            for p in g["parts"]:
                if p[0] == "doctest":
                    last_py = norm(p[2])
        elif last_py is not None:
            parts = []
            for p in g["parts"]:
                if p[0] == "text":
                    parts += [["text", t] for t in "\n".join(p[1]).strip().split("\n\n") if t.strip()]
                else:
                    parts.append([p[0], p[-1].strip("\n")])
            out.setdefault(last_py, []).append([g["lang"], parts])
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
