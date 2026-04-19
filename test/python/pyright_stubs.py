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
      summary = json.loads(result.stdout)["summary"]
    except (ValueError, KeyError) as e:
      self.fail(
          "could not parse pyright --outputjson: %s\nstdout:\n%s\nstderr:\n%s"
          % (e, result.stdout[:2000], result.stderr[:2000])
      )
    errors = summary["errorCount"]
    if errors > PYRIGHT_ERROR_BUDGET:
      # Surface the first few diagnostics to make regressions actionable.
      diagnostics = json.loads(result.stdout).get("generalDiagnostics", [])
      sample = "\n".join(
          "  %s:%s  %s" % (d.get("file", "?"),
                           d.get("range", {}).get("start", {}).get("line", "?"),
                           d.get("message", "").splitlines()[0])
          for d in diagnostics if d.get("severity") == "error"
      )
      self.fail(
          "pyright errors regressed: %d > budget %d\nFirst errors:\n%s"
          % (errors, PYRIGHT_ERROR_BUDGET, sample[:4000])
      )
    if errors < PYRIGHT_ERROR_BUDGET:
      print(
          "pyright error count improved: %d < budget %d -- "
          "consider lowering PYRIGHT_ERROR_BUDGET in %s"
          % (errors, PYRIGHT_ERROR_BUDGET, os.path.basename(__file__)),
          file=sys.stderr,
      )


if __name__ == "__main__":
  unittest.main()
