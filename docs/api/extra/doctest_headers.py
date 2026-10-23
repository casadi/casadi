#!/usr/bin/env python
"""Run the \\doctest blocks of CasADi C++ headers directly, without regenerating doc.i.

Usage: python docs/api/extra/doctest_headers.py [--verbose] [--lang python|matlab|julia|javascript|all] [header.hpp ...]
With no headers all headers under casadi/ are scanned. Default language: python.
Block syntax in a header comment:
  \\doctest             Python, doctest format (>>> lines + expected output)
  \\doctest{julia}      Documenter format (julia> lines + expected output), run by doctest_julia.jl
  \\doctest{javascript} node REPL format (> lines + expected output), run by doctest_js.js
  \\doctest{matlab}     command-window format (>> lines + expected output), run by doctest_matlab.m
                       in MATLAB (install-matlab/matlab or CASADI_MATLAB_PATH; MATLAB_CMD=octave for Octave)
  \\enddoctest
Needs an importable casadi on PYTHONPATH (python), build/swig/julia/CasADiNative or CASADI_JL (julia),
build-wasm/swig/wasm-js/casadi.js or CASADI_WASM_JS (javascript).
"""
import doctest
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
JULIA = os.environ.get("JULIA", "julia")


def collect(path):
    """Yield (lineno, lang, text) for every \\doctest ... \\enddoctest block in a header."""
    blocks = []
    with open(path) as f:
        lines = f.read().split("\n")
    cur = None
    for i, line in enumerate(lines, 1):
        s = re.sub(r"^\s*\*?\s?", "", line)
        m = re.match(r"\\doctest(?:\{(\w+)\})?\s*$", s.strip())
        if m:
            if cur is not None:
                raise ValueError("%s:%d: nested \\doctest" % (path, i))
            cur = (i, m.group(1) or "python", [])
        elif "\\enddoctest" in s:
            if cur is None:
                raise ValueError("%s:%d: \\enddoctest without \\doctest" % (path, i))
            blocks.append((cur[0], cur[1], "\n".join(cur[2]) + "\n"))
            cur = None
        elif "\\doctestout" in s or "\\doctest" in s:
            raise ValueError("%s:%d: malformed \\doctest line: %r" % (path, i, line))
        elif cur is not None:
            cur[2].append(s)
    if cur is not None:
        raise ValueError("%s:%d: unterminated \\doctest" % (path, cur[0]))
    return blocks


def make_globs():
    import casadi
    import numpy as np
    globs = {"__name__": "casadi_doctest", "casadi": casadi, "ca": casadi, "np": np, "numpy": np}
    exec("from casadi import *", globs)
    return globs


def run_python(blocks, verbose):
    parser = doctest.DocTestParser()
    flags = doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE
    runner = doctest.DocTestRunner(verbose=verbose, optionflags=flags)
    for name, text in blocks:
        path, lineno = name.rsplit(":", 1)
        runner.run(parser.get_doctest(text, make_globs(), name, path, int(lineno)))
    failed, attempted = runner.summarize(verbose=False)
    print("python: %d blocks, %d examples, %d failures" % (len(blocks), attempted, failed))
    return failed == 0


def run_external(lang, blocks, verbose):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        for name, text in blocks:
            f.write("#### %s\n%s\n" % (name, text))
        tmp = f.name
    if lang == "julia":
        cmd = [JULIA, os.path.join(HERE, "doctest_julia.jl"), tmp]
    else:
        cmd = ["node", os.path.join(HERE, "doctest_js.js"), tmp]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out = p.stdout.strip()
    print(out if verbose or p.returncode else out.split("\n")[-1])
    os.unlink(tmp)
    return p.returncode == 0


def run_matlab(files, verbose):
    """Run doctest_matlab.m in MATLAB (default) or Octave (MATLAB_CMD=octave)."""
    path = os.environ.get("CASADI_MATLAB_PATH", os.path.join(ROOT, "install-matlab", "matlab"))
    call = "addpath('%s'); addpath('%s'); n = doctest_matlab(%s%s); exit(double(n > 0))" % (
        path, HERE, ", ".join("'%s'" % f for f in files), ", '-verbose'" if verbose else "")
    env = dict(os.environ)
    if os.environ.get("MATLAB_CMD", "matlab") == "octave":
        cmd = ["octave", "--no-gui", "--no-window-system", "--eval", call]
    else:
        cmd = ["matlab", "-nodisplay", "-nosplash", "-nodesktop", "-batch", call]
        if os.path.exists("/usr/lib/x86_64-linux-gnu/libstdc++.so.6"):  # MATLAB bundles an older one
            env["LD_PRELOAD"] = "/usr/lib/x86_64-linux-gnu/libstdc++.so.6"
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
    out = "\n".join(l for l in p.stdout.split("\n") if "import has the same name" not in l
                    and not l.startswith(("[> In ", "In addpath", "In path", "warning: load_path")))
    print(out.strip() if verbose or p.returncode else out.strip().split("\n")[-1])
    return p.returncode == 0


def run(paths, langs, verbose=False):
    by_lang = {}
    for path in paths:
        for lineno, lang, text in collect(path):
            name = os.path.relpath(path, ROOT) if os.path.abspath(path).startswith(ROOT) else path
            by_lang.setdefault(lang, []).append(("%s:%d" % (name, lineno), text))
    ok = True
    for lang in sorted(by_lang):
        if "all" not in langs and lang not in langs:
            print("%s: %d blocks skipped (use --lang %s)" % (lang, len(by_lang[lang]), lang))
            continue
        if lang == "python":
            ok &= run_python(by_lang[lang], verbose)
        elif lang in ("julia", "javascript"):
            ok &= run_external(lang, by_lang[lang], verbose)
        elif lang == "matlab":
            files = sorted({name.rsplit(":", 1)[0] for name, _ in by_lang[lang]})
            ok &= run_matlab([f if os.path.isabs(f) else os.path.join(ROOT, f) for f in files], verbose)
        else:
            raise ValueError("unknown doctest language %r" % lang)
    if not by_lang:
        print("doctest_headers: no blocks found")
    return ok


if __name__ == "__main__":
    args = sys.argv[1:]
    verbose = "--verbose" in args
    langs = ["python"]
    if "--lang" in args:
        i = args.index("--lang")
        langs = args[i + 1].split(",")
        del args[i:i + 2]
    args = [a for a in args if not a.startswith("--")]
    if not args:
        for dirpath, dirnames, files in os.walk(os.path.join(ROOT, "casadi")):
            dirnames[:] = [d for d in dirnames if d != "include"]  # installed copies of the headers
            args += [os.path.join(dirpath, f) for f in files if f.endswith(".hpp")]
    sys.exit(0 if run(sorted(args), langs, verbose) else 1)
