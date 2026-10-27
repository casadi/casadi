#!/usr/bin/env python
"""Run the \\doctest blocks of CasADi C++ headers directly, without regenerating doc.i.

Usage: python docs/api/extra/doctest_headers.py [--verbose] [--lang python|matlab|julia|javascript|java|all] [header.hpp ...]
With no headers all headers under casadi/ are scanned. Default language: python.
Block syntax in a header comment:
  \\lang{julia}          everything up to \\endlang (prose, lists, examples) is for this language only
  \\doctest              example in the format of the enclosing \\lang group, which it requires:
                        python: doctest format (>>> lines + expected output)
                        julia: Documenter format (julia> lines + expected output), run by doctest_julia.jl
                        javascript: node REPL format (> lines + expected output), run by doctest_js.js
                        matlab: command-window format (>> lines + expected output), run by doctest_matlab.m
                        in MATLAB (install-matlab/matlab or CASADI_MATLAB_PATH; MATLAB_CMD=octave for Octave)
                        java: jshell format (jshell> lines, "...> " continuations, expected output
                        without jshell's "$1 ==> " prefix), run in jshell with org.casadi.* and java.util.*
                        imported; CASADI_JAVA_PATH is the directory of casadi.jar (default install/lib)
  \\enddoctest
  \\endlang
Numbers: an expected number matches the actual one at the same place when they differ by at most
half a unit in its last written digit (0 for an integer) or 1e-9 relative (absolute for 0),
whichever is larger: 0.3679 matches 0.367879, 0 matches 1e-17. Applies when both outputs hold
equally many numbers; 00 (structural zero) is text. Same rule in every runner.
Needs an importable casadi on PYTHONPATH (python), an importable CasADiNative (julia; see
test/julia/package/make_package.jl),
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


def parse(path):
    """Return the \\lang groups of a header: dicts with lang, line, comment and parts, a list of
    ("text", lines), ("verbatim", text) and ("doctest", lineno, text) in source order."""
    with open(path) as f:
        lines = f.read().split("\n")
    groups, group, block, comment = [], None, None, 0
    for i, line in enumerate(lines, 1):
        if "/**" in line or "/*!" in line:
            comment += 1
        s = re.sub(r"^\s*\*?\s?", "", line)
        t = re.sub(r"\*/\s*$", "", s).strip()  # '\enddoctest */' closes the comment too
        where = "%s:%d" % (path, i)
        m = re.match(r"\\lang\{(\w+)\}$", t)
        if m:
            if group is not None:
                raise ValueError("%s: nested \\lang group" % where)
            group = {"lang": m.group(1), "line": i, "comment": comment, "parts": []}
        elif t == "\\endlang":
            if group is None or block is not None:
                raise ValueError("%s: \\endlang without \\lang" % where)
            groups.append(group)
            group = None
        elif t == "\\doctest" or (t == "\\verbatim" and group is not None):
            if group is None:
                raise ValueError("%s: \\doctest outside a \\lang group" % where)
            if block is not None:
                raise ValueError("%s: nested %s" % (where, t))
            block = (t[1:], i, [])
            indent = line[:line.index(t)]  # julia continuation lines keep their 7 spaces
        elif block is not None and t == "\\end" + block[0]:
            text = "\n".join(block[2]) + "\n"
            group["parts"].append(("doctest", block[1], text) if block[0] == "doctest" else ("verbatim", text))
            block = None
        elif "\\lang" in s or "\\endlang" in s or "\\doctest" in s:
            raise ValueError("%s: malformed line: %r" % (where, line))
        elif block is not None:
            block[2].append(line[len(indent):] if group["lang"] == "julia" and line.startswith(indent) else s)
        elif group is not None:
            if group["parts"] and group["parts"][-1][0] == "text":
                group["parts"][-1][1].append(t)
            else:
                group["parts"].append(("text", [t]))
        if "*/" in line and (group is not None or block is not None):
            raise ValueError("%s: unterminated \\%s" % (where, "lang" if block is None else block[0]))
    return groups


def collect(path):
    """Return (lineno, lang, text) for every \\doctest ... \\enddoctest block in a header."""
    return [(p[1], g["lang"], p[2]) for g in parse(path) for p in g["parts"] if p[0] == "doctest"]


NUMBER = re.compile(r"(?<![\w.@])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?(?![\w.])")


def numbers(s):
    return [m for m in NUMBER.finditer(s) if not re.match(r"[-+]?0\d", m.group())]


def tolerance(lit):
    m = re.match(r"[-+]?\d*(?:\.(\d*))?(?:[eE]([-+]?\d+))?$", lit)
    tol = 0.5 * 10.0 ** (int(m.group(2) or 0) - len(m.group(1) or "")) if m.group(1) is not None or m.group(2) else 0.0
    return max(tol, 1e-9 * (abs(float(lit)) or 1.0))


def as_expected(want, got):
    """got with each number that matches its expected counterpart written as expected"""
    w, g = numbers(want), numbers(got)
    if len(w) != len(g):
        return got
    out, pos = [], 0
    for a, b in zip(w, g):
        if abs(float(a.group()) - float(b.group())) <= tolerance(a.group()) * (1 + 1e-9):
            out += [got[pos:b.start()], a.group()]
            pos = b.end()
    return "".join(out) + got[pos:]


class Checker(doctest.OutputChecker):
    def check_output(self, want, got, optionflags):
        return super().check_output(want, as_expected(want, got), optionflags)


def make_globs():
    import casadi
    import numpy as np
    globs = {"__name__": "casadi_doctest", "casadi": casadi, "ca": casadi, "np": np, "numpy": np}
    exec("from casadi import *", globs)
    return globs


def run_python(blocks, verbose):
    parser = doctest.DocTestParser()
    flags = doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE
    runner = doctest.DocTestRunner(Checker(), verbose=verbose, optionflags=flags)
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


JAVA_MARK = "@@casadi-doctest@@"


def java_examples(text):
    """(source, expected) pairs of a jshell-format block"""
    out = []
    for line in text.split("\n"):
        if line.startswith("jshell> "):
            out.append([line[8:], []])
        elif re.match(r"\s*\.\.\.> ", line) and out and not out[-1][1]:
            out[-1][0] += "\n" + re.sub(r"^\s*\.\.\.> ", "", line)
        elif line.strip() and out:
            out[-1][1].append(line)
    return [(src, "\n".join(want)) for src, want in out]


def run_java(blocks, verbose):
    """Run all blocks in one jshell session, /reset between blocks"""
    path = os.environ.get("CASADI_JAVA_PATH", os.path.join(ROOT, "install", "lib"))
    mark = 'System.out.println("%s")' % JAVA_MARK
    script, examples = [], []
    for name, text in blocks:
        ex = java_examples(text)
        examples.append((name, ex))
        # jshell's startup imports java.util.function.*: org.casadi.Function wins
        script += ["/reset", "import org.casadi.*;", "import org.casadi.Function;", "import java.util.*;", mark]
        for src, _ in ex:
            script += [src, mark]
    script.append("/exit")
    cmd = ["jshell", "-q", "-R-Dcasadi.library.path=" + path, "--class-path", os.path.join(path, "casadi.jar")]
    # stderr (casadi warnings) is not compared, as in Python's doctest
    p = subprocess.run(cmd, input="\n".join(script) + "\n", stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    chunks = [re.sub(r"(?m)^\$\d+ ==> ", "", re.sub(r"^[A-Za-z_]\w* ==> ", "", c.replace("jshell> ", "").strip("\n")))
              for c in p.stdout.split(JAVA_MARK + "\n")]
    checker = Checker()
    flags = doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE
    failed, attempted, k = 0, 0, 0
    for name, ex in examples:
        k += 1  # output of the block preamble
        for src, want in ex:
            got = chunks[k] if k < len(chunks) else "<jshell output ended>"
            k += 1
            attempted += 1
            if not checker.check_output(want + "\n", got + "\n", flags):
                failed += 1
                print("%s\nFailed example:\n    %s\nExpected:\n    %s\nGot:\n    %s\n" % (
                    name, src.replace("\n", "\n    "), want.replace("\n", "\n    "), got.replace("\n", "\n    ")))
            elif verbose:
                print("%s: ok: %s" % (name, src.split("\n")[0]))
    print("java: %d blocks, %d examples, %d failures" % (len(blocks), attempted, failed))
    return failed == 0


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
        elif lang == "java":
            ok &= run_java(by_lang[lang], verbose)
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
