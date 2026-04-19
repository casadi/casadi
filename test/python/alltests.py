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
import importlib
import sys
import unittest

# Each submodule owns its TestCase subclasses; we import them as
# namespaces and let TestLoader discover their tests.  Avoids
# `from X import *` polluting the global scope with every symbol from
# every test module (which also caused pyright Final-reassignment
# cascades on names like `inf`/`pi`).
_TEST_MODULES = [
    "mx", "sx", "typemaps", "integration", "ocp", "nlp",
    "implicitfunction", "ad", "sparsity", "linearsolver", "matrix",
    "conic", "misc", "function", "tools", "simulator", "vectortools",
    "optistack", "feasiblesqpmethod", "serialize", "threads",
    "daebuilder", "pyright_stubs",
]


def build_suite():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    for name in _TEST_MODULES:
        mod = importlib.import_module(name)
        suite.addTests(loader.loadTestsFromModule(mod))
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(build_suite())
    sys.exit(0 if result.wasSuccessful() else 1)
