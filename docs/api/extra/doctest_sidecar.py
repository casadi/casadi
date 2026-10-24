#!/usr/bin/env python
"""Pair the non-Python \\lang groups of the C++ headers with the Python example they translate.

Writes JSON {normalized python block: [[lang, parts], ...]} for the Sphinx hook, which inserts
those groups after the matching Python example on the API page. parts is a list of
["text", paragraph text], ["verbatim", text] and ["doctest", text]. A group in another language
belongs to the closest preceding Python group of the same Doxygen comment and follows its last
\\doctest; groups without such an anchor are left out. A first ["python", parts] entry lists the
paragraphs of the Python group itself before its first \\doctest, which only apply to Python.
Usage: python doctest_sidecar.py <casadi source dir> <output.json>
"""
import json
import os
import sys

from doctest_headers import parse


def norm(text):
    return "\n".join(l.strip() for l in text.strip().split("\n") if l.strip())


def paragraphs(lines):
    return [["text", t] for t in "\n".join(lines).strip().split("\n\n") if t.strip()]


def pairs(path):
    out = {}
    comment, last_py = None, None
    for g in parse(path):
        if g["comment"] != comment:
            comment, last_py = g["comment"], None
        if g["lang"] == "python":
            last_py, py_text = None, []
            for p in g["parts"]:
                if p[0] == "doctest":
                    last_py = norm(p[2])
                elif p[0] == "text" and last_py is None:
                    py_text += paragraphs(p[1])
        elif last_py is not None:
            parts = []
            for p in g["parts"]:
                if p[0] == "text":
                    parts += paragraphs(p[1])
                else:
                    parts.append([p[0], p[-1].strip("\n")])
            entry = out.setdefault(last_py, [])
            if py_text and not entry:
                entry.append(["python", py_text])
            entry.append([g["lang"], parts])
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
