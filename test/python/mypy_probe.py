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
"""Type-checked against the generated stubs by CI; never executed as a test.

mypy silences everything *inside* a PEP 561 package, so this file cannot
police stub quality in general -- pyright_stubs.py does that.  What it does
catch is the class of defect that is NOT silenceable: anything that stops
mypy building casadi.pyi's AST is a blocking error which no user config can
ignore, and it takes the whole package down with it (issue #4390 shipped a
duplicate `INOUT` parameter and forced users to pin casadi<3.8).

So: exercise enough of the surface that mypy has to load the stubs, and keep
one call per stub construct that has been broken before.
"""
from typing import TYPE_CHECKING

# Everything below is for the type checker only.  This directory's runners
# (test/internal/testsuite.py) execute every .py they find, and some of these
# calls are not safe to actually run -- Sparsity.dfs aborts on out-of-contract
# vector sizes.  mypy and pyright both analyse TYPE_CHECKING blocks, so the
# file stays fully checked while being a no-op when executed.  Everything the
# checks need is imported in here too, so no cell has to ship typing_extensions
# just to walk this directory.
if TYPE_CHECKING:
  import numpy as np
  from typing_extensions import assert_type

  import casadi as ca

  # Recursive options dict: bool/int/float/str/Function/None at any depth.
  x = ca.SX.sym("x")
  f = ca.Function("f", [x], [x ** 2])
  solver = ca.nlpsol("solver", "ipopt", {"x": x, "f": x * x},
                     {"ipopt": {"tol": 1e-6}, "print_time": True,
                      "iteration_callback": f, "unused": None})
  stats = solver.stats()

  # Options round-trip, including the None entries Options.sanitize preserves.
  sanitized = ca.Options.sanitize({"foo": 1, "bar": {"r": None}})

  # numpy scalars are accepted by the DM operators.
  dm = ca.DM([[1, 2], [3, 4]])
  scaled = np.float64(2.0) - dm

  # &INOUT typemaps: parameter names must be distinct, return types must be types.
  sp = ca.Sparsity.diag(3)
  nz = sp.get_nz([0, 1, 2])
  dedup = sp.removeDuplicates([0, 1, 2])
  top, xi, pstack, marked = sp.dfs(0, 0, [0, 0, 0], [0, 0, 0], [0, 1, 2],
                                   [False, False, False])
  subs = ca.substitute_inplace([x], [x], [x], False)

  # IM-typed return (MX::mapping()) resolves to DM.
  mx = ca.MX.sym("mx", 2, 2)
  mapping = mx[0, 0].mapping()

  # assert_type pins the resolved types.  Importing casadi is not enough to
  # police the stubs: mypy suppresses diagnostics inside a PEP 561 package,
  # so a construct whose type silently degrades produces no error at all at
  # the call site.  These assertions fail instead.
  assert_type(mapping, ca.DM)                      # not the undefined name IM
  assert_type(nz, list[int])                       # not the doc token [int]
  assert_type(dedup, list[int])
  assert_type(mx[mx], ca.MX)                       # MX-valued index (MX_get)
  assert_type(mx[0, 0], ca.MX)
  assert_type(dm[0, 0], ca.DM)
  assert_type(ca.SX.sym("s", 2, 2)[0], ca.SX)
  assert_type(scaled, ca.DM)                       # numpy scalar operand
  opti = ca.Opti()
  assert_type(opti.variable(2, 2, "symmetric"), ca.MX)   # builtins.str arg
  mx[0, 0] = 1
  mx[mx] = 1
