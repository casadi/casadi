#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             KU Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
from __future__ import print_function  # `print(..., file=...)` below requires this under Py2

import json
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest

import casadi
from helpers import *

_PY310 = sys.version_info >= (3, 10)


# Error budget for the full-suite pyright run (test_pyright_suite).
# Ratchet DOWN as stub quality improves; never raise it silently.
PYRIGHT_ERROR_BUDGET = 0


def _have_pyright():
  return shutil.which("pyright") is not None


def _casadi_package_dir():
  return os.path.dirname(os.path.abspath(casadi.__file__))


import re
import tokenize


# Match `# expect-error` optionally followed by `: rule1, rule2, ...`.
# Rule names are reportXxxxx identifiers.  Anchoring the rule list to
# start with an identifier character rejects placeholder text like
# `# expect-error: <rule1>,<rule2>` (the angle brackets aren't valid).
_EXPECT_ERROR_RE = re.compile(
    r"^\s*expect-error\b(?:\s*:\s*([A-Za-z_][A-Za-z_0-9, ]*))?\s*$"
)


def _scan_expect_errors(tracked_files):
  """Walk tracked .py files; return {(realpath, line_1based): set(rule)}.

  A sentinel comment `# expect-error: <rule1>,<rule2>` on a Python
  line asserts that pyright SHOULD report at least one diagnostic on
  that line whose `rule` is one of the listed rules.  An empty rule
  list (`# expect-error`) accepts any rule.

  The sentinel converts what would otherwise be a `# pyright: ignore`
  suppression into an explicit positive test: if the stubs become too
  permissive and the previously-rejected expression now type-checks,
  the test fails instead of silently passing.

  Uses Python's tokenize module so the scanner only sees COMMENT
  tokens, not text inside docstrings / regex source strings that
  happen to mention the sentinel syntax.
  """
  sentinels = {}
  for f in tracked_files:
    if not f.endswith(".py"):
      continue
    try:
      with open(f, "rb") as fh:
        tokens = list(tokenize.tokenize(fh.readline))
    except (OSError, tokenize.TokenError, SyntaxError):
      continue
    for tok in tokens:
      if tok.type != tokenize.COMMENT:
        continue
      # tok.string includes the leading `#`; strip it before matching.
      body = tok.string.lstrip("#")
      m = _EXPECT_ERROR_RE.match(body)
      if m is None:
        continue
      raw = (m.group(1) or "").strip()
      rules = set()
      for r in raw.split(","):
        r = r.strip()
        if re.match(r"^[A-Za-z_][A-Za-z_0-9]*$", r):
          rules.add(r)
      try:
        key = (os.path.realpath(f), tok.start[0])
      except OSError:
        continue
      sentinels[key] = rules
  return sentinels


def _index_error_diagnostics(diagnostics):
  """Group error-severity pyright diagnostics by (realpath, line_1based)
  -> list of raw diagnostic dicts.  Pyright line numbers are 0-based;
  we normalise to 1-based to match grep / editor conventions."""
  idx = {}
  for d in diagnostics:
    if d.get("severity") != "error":
      continue
    f = d.get("file") or ""
    line0 = d.get("range", {}).get("start", {}).get("line")
    if not f or line0 is None:
      continue
    try:
      key = (os.path.realpath(f), int(line0) + 1)
    except OSError:
      continue
    idx.setdefault(key, []).append(d)
  return idx


def _tracked_files_in(directory):
  """Set of realpath()'d files tracked by git under `directory`.
  Returns None if git is unavailable or `directory` is not inside a
  git working tree (CI clean checkouts, source-tarball builds)."""
  try:
    root = subprocess.check_output(
        ["git", "rev-parse", "--show-toplevel"],
        cwd=directory, text=True, stderr=subprocess.DEVNULL,
    ).strip()
    listing = subprocess.check_output(
        ["git", "ls-files", "--full-name", "--", directory],
        cwd=root, text=True, stderr=subprocess.DEVNULL,
    )
  except (OSError, subprocess.CalledProcessError):
    return None
  out = set()
  for rel in listing.splitlines():
    rel = rel.strip()
    if not rel:
      continue
    try:
      out.add(os.path.realpath(os.path.join(root, rel)))
    except OSError:
      pass
  return out


SMOKE_SRC = textwrap.dedent(
    """
    # assert_type arrived in typing in Py3.11; typing_extensions
    # backports it.  Using that unconditionally keeps the smoke test
    # portable across all supported Python targets.
    from typing_extensions import assert_type

    import casadi as ca

    x = ca.MX.sym("x")
    y = ca.sin(x)
    assert_type(x, ca.MX)
    assert_type(y, ca.MX)
    assert_type(ca.vertcat(x, y), ca.MX)

    sx = ca.SX.sym("sx")
    assert_type(ca.gradient(ca.dot(sx, sx), sx), ca.SX)

    dense = ca.DM.eye(3)
    assert_type(dense, ca.DM)
    assert_type(ca.mtimes(dense, dense), ca.DM)

    f = ca.Function("f", [x], [y])
    assert_type(f, ca.Function)

    solver = ca.nlpsol("solver", "ipopt", {"x": x, "f": y}, {})
    assert_type(solver, ca.Function)

    opti = ca.Opti()
    decision = opti.variable()
    assert_type(decision, ca.MX)
    assert_type(opti.solve(), ca.OptiSol)
    """
)

# Expressions that MUST be rejected by the stubs.  Each entry is a
# single-line Python snippet.  The negative-smoke test wraps each in
# `reveal_type(EXPR)` and asserts pyright reports at least one error on
# that line.  Used to guard against stub regressions that silently
# accept meaningless casadi calls.
NEGATIVE_CASES = [
    # Mixing symbolic types at the operator level is a runtime error in
    # casadi (to_ptr<MX>(SX) / to_ptr<SX>(MX) fail); stubs must reject it.
    "ca.MX.sym('x') + ca.SX.sym('y')",
    "ca.SX.sym('x') + ca.MX.sym('y')",
    "ca.MX.sym('x') * ca.SX.sym('y')",
    "ca.MX.sym('x') - ca.SX.sym('y')",
]


@unittest.skipUnless(_PY310, "pyright_stubs requires Python 3.10+ (pyright itself does not support Python 2)")
class TypingTests(casadiTestCase):

  def test_stubs_installed(self):
    """Generated stub files are present next to the package."""
    pkg = _casadi_package_dir()
    for f in ("casadi.pyi", "__init__.pyi", "py.typed"):
      path = os.path.join(pkg, f)
      self.assertTrue(
          os.path.exists(path),
          "expected stub file %s missing; is -stubs enabled in the SWIG invocation?"
          % path,
      )

  def test_pyright_smoke(self):
    """pyright accepts a curated smoke script against the installed stubs."""
    if not _have_pyright():
      self.skipTest("pyright not installed")
    with tempfile.TemporaryDirectory() as td:
      src = os.path.join(td, "smoke.py")
      with open(src, "w") as fh:
        fh.write(SMOKE_SRC)
      # Point pyright at the installed casadi package so it uses our stubs.
      env = os.environ.copy()
      env.setdefault("PYTHONPATH", os.path.dirname(_casadi_package_dir()))
      # Let pyright use the running interpreter's version (so the same
      # test suite exercises Py3.9/3.10/.../3.13 stubs end-to-end).
      result = subprocess.run(
          ["pyright", src],
          capture_output=True, text=True, env=env,
      )
      self.assertEqual(
          result.returncode, 0,
          "pyright smoke failed:\nstdout:\n%s\nstderr:\n%s"
          % (result.stdout, result.stderr),
      )

  def test_pyright_negative(self):
    """Expressions in NEGATIVE_CASES must each be rejected by pyright.

    Guards against stubs silently accepting known-bad patterns (e.g.
    mixing MX with SX) after a future stub refactor widens a union too
    far.  The test builds one file with each case on its own line,
    then asserts pyright reports an error on exactly those lines.
    """
    if not _have_pyright():
      self.skipTest("pyright not installed")
    with tempfile.TemporaryDirectory() as td:
      src = os.path.join(td, "negatives.py")
      lines = ["import casadi as ca"]
      for case in NEGATIVE_CASES:
        lines.append("_ = " + case)
      with open(src, "w") as fh:
        fh.write("\n".join(lines) + "\n")
      env = os.environ.copy()
      env.setdefault("PYTHONPATH", os.path.dirname(_casadi_package_dir()))
      result = subprocess.run(
          ["pyright", "--outputjson", src],
          capture_output=True, text=True, env=env,
      )
      try:
        diagnostics = json.loads(result.stdout).get("generalDiagnostics", [])
      except ValueError as e:
        self.fail("could not parse pyright output: %s\n%s" % (e, result.stdout[:2000]))
      # Map error-bearing lines.  Line numbers are 0-based in pyright's
      # JSON; our cases start on line index 1 (after the import).
      error_lines = {
          d["range"]["start"]["line"]
          for d in diagnostics if d.get("severity") == "error"
      }
      for idx, case in enumerate(NEGATIVE_CASES, start=1):
        self.assertIn(
            idx, error_lines,
            "stubs silently accept '%s' -- expected pyright to reject it.\n"
            "Full pyright stdout:\n%s" % (case, result.stdout[:2000]),
        )

  @memory_heavy()
  def test_pyright_suite(self):
    """pyright on the full python test suite stays within the error budget."""
    if not _have_pyright():
      self.skipTest("pyright not installed")
    if not _PY310:
      self.skipTest("PYRIGHT_ERROR_BUDGET is calibrated against Py3.10+ typeshed")
    env = os.environ.copy()
    env.setdefault("PYTHONPATH", os.path.dirname(_casadi_package_dir()))
    test_dir = os.path.dirname(os.path.abspath(__file__))
    # Run from the test dir so sibling imports ("from helpers import *",
    # "import mx", etc.) resolve; otherwise pyright reports ~700 spurious
    # reportMissingImports/reportUndefinedVariable diagnostics.
    result = subprocess.run(
        ["pyright", "--outputjson", "."],
        capture_output=True, text=True, env=env, cwd=test_dir,
    )
    try:
      diagnostics = json.loads(result.stdout).get("generalDiagnostics", [])
    except (ValueError, KeyError) as e:
      self.fail(
          "could not parse pyright --outputjson: %s\nstdout:\n%s\nstderr:\n%s"
          % (e, result.stdout[:2000], result.stderr[:2000])
      )

    # Count errors ONLY in git-tracked files.  Locally-developers' clutter
    # under test/python/ (debug_*.py scripts, scratch numpy_experiments
    # variants, etc.) would otherwise bust the budget for no good reason
    # -- CI never sees those files because they're not committed.
    tracked = _tracked_files_in(test_dir)
    def _in_tracked(d):
      if tracked is None:
        return True   # couldn't enumerate; be inclusive
      f = d.get("file") or ""
      try:
        return os.path.realpath(f) in tracked
      except OSError:
        return False
    # `# expect-error: <rule>` sentinels mark intentional negative type
    # tests: the diagnostic on that line is expected.  Drop those from
    # the budget count; the separate test_pyright_expect_error_sentinels
    # test verifies they actually fired.
    sentinel_keys = set(
        _scan_expect_errors(tracked).keys() if tracked else ()
    )
    diag_idx = _index_error_diagnostics(diagnostics)
    error_diagnostics = []
    for d in diagnostics:
      if d.get("severity") != "error":
        continue
      if not _in_tracked(d):
        continue
      f = d.get("file") or ""
      line0 = d.get("range", {}).get("start", {}).get("line")
      if line0 is not None:
        try:
          if (os.path.realpath(f), int(line0) + 1) in sentinel_keys:
            continue
        except OSError:
          pass
      error_diagnostics.append(d)
    errors = len(error_diagnostics)
    if errors > PYRIGHT_ERROR_BUDGET:
      # Surface the first few diagnostics to make regressions actionable.
      sample = "\n".join(
          "  %s:%s  %s" % (d.get("file", "?"),
                           d.get("range", {}).get("start", {}).get("line", "?"),
                           d.get("message", "").splitlines()[0])
          for d in error_diagnostics
      )
      self.fail(
          "pyright errors regressed: %d > budget %d (tracked files only)\n"
          "First errors:\n%s"
          % (errors, PYRIGHT_ERROR_BUDGET, sample[:4000])
      )
    if errors < PYRIGHT_ERROR_BUDGET:
      print(
          "pyright error count improved: %d < budget %d -- "
          "consider lowering PYRIGHT_ERROR_BUDGET in %s"
          % (errors, PYRIGHT_ERROR_BUDGET, os.path.basename(__file__)),
          file=sys.stderr,
      )

  @memory_heavy()
  def test_pyright_expect_error_sentinels(self):
    """Each `# expect-error: <rule1>,<rule2>` sentinel in a tracked
    file must correspond to a pyright diagnostic on that line whose
    rule matches one of the listed rules.

    The sentinel is the type-system analogue of `assertRaises` -- it
    pins that a particular expression MUST fail type-checking.  If a
    later stub widening silently accepts the expression, this test
    fails before the change ships, instead of the runtime `assertRaises`
    quietly passing while users get type-clean rope to hang themselves
    with.
    """
    if not _have_pyright():
      self.skipTest("pyright not installed")
    if not _PY310:
      self.skipTest("PYRIGHT_ERROR_BUDGET is calibrated against Py3.10+ typeshed")
    env = os.environ.copy()
    env.setdefault("PYTHONPATH", os.path.dirname(_casadi_package_dir()))
    test_dir = os.path.dirname(os.path.abspath(__file__))
    tracked = _tracked_files_in(test_dir)
    if not tracked:
      self.skipTest("not in a git working tree; cannot enumerate "
                    "tracked files for sentinel scan")
    sentinels = _scan_expect_errors(tracked)
    if not sentinels:
      # No sentinels in the codebase yet; nothing to verify.
      return
    result = subprocess.run(
        ["pyright", "--outputjson", "."],
        capture_output=True, text=True, env=env, cwd=test_dir,
    )
    try:
      diagnostics = json.loads(result.stdout).get("generalDiagnostics", [])
    except ValueError as e:
      self.fail("could not parse pyright --outputjson: %s\n%s"
                % (e, result.stdout[:2000]))
    diag_idx = _index_error_diagnostics(diagnostics)
    failures = []
    for (f, line), expected_rules in sorted(sentinels.items()):
      diags = diag_idx.get((f, line), [])
      if not diags:
        failures.append(
            "  %s:%d -- expected pyright error %s but found none. "
            "Stubs may have become too permissive."
            % (f, line, sorted(expected_rules) or "<any rule>")
        )
        continue
      if expected_rules:
        actual = {d.get("rule") for d in diags if d.get("rule")}
        if not (expected_rules & actual):
          failures.append(
              "  %s:%d -- expected rule(s) %s but got %s"
              % (f, line, sorted(expected_rules),
                 sorted(r for r in actual))
          )
    if failures:
      self.fail(
          "expect-error sentinels not satisfied (%d failures):\n%s"
          % (len(failures), "\n".join(failures[:50]))
      )


if __name__ == "__main__":
  unittest.main()
