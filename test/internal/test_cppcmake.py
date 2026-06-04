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
import sys
src = sys.argv[1]
import os

"""
Test cpp files that are made with CMAKE and put in the /build/bin directory.

By default the executables are run from ../build/bin. Set CASADI_EXAMPLE_BINDIR
to run pre-installed executables from another directory (e.g. an unpacked
binary package); sources without a matching executable there are skipped.
"""

from testsuite import TestSuite

from subprocess import *

args = sys.argv[2:]
bindir = os.environ.get("CASADI_EXAMPLE_BINDIR")
if bindir:
  workingdir = lambda dir: bindir
  command = lambda dir,fn,opt: ['./'+fn.replace('.cpp','')]+opt
  missing = [f for f in os.listdir(src) if f.endswith('.cpp')
             and not any(os.path.exists(os.path.join(bindir, f[:-4]+ext)) for ext in ('', '.exe'))]
  if missing:
    args = args + ['-skipfiles='+' '.join(missing)]
else:
  workingdir = lambda dir: '../build'
  command = lambda dir,fn,opt: ['./bin/'+fn.replace('.cpp','')]+opt

t = TestSuite(dirname=src,
  suffix="cpp",
  command = command,
  workingdir = workingdir,
  skipdirs=[".svn","ctemplate","cmake_find_package","cmake_pkgconfig","CMakeFiles"],
  args=args
  )

t.run()




