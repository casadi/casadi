#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
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
import shutil

from testsuite import TestSuite

def setdummybackend(dir):
  """
  Set the matplotlib backend to a dummy value, so that no graphs are generated.
  """
  #shutil.copy('internal/pylab.py',dir)
  f=open(dir+'/matplotlibrc','w')
  f.write('backend : Template')
  f.close()

def removedummybackend(dir):
  if os.path.exists(dir+'/matplotlibrc'):
    os.remove(dir+'/matplotlibrc')
  if os.path.exists(dir+'/pylab.pyc'):
    os.remove(dir+'/pylab.pyc')
  if os.path.exists(dir+'/pylab.py'):
    os.remove(dir+'/pylab.py')

python = "python"

if "WITH_PYTHON3" in os.environ:
  python = "python3"

t = TestSuite(dirname=src,
  suffix="py",
  preRun=setdummybackend,
  postRun=removedummybackend,
  command = lambda dir,fn, opt:  [python,"-W","error::SyntaxWarning","-W","error:This CasADi:DeprecationWarning", fn] + opt,
  skipdirs=[".svn","ctemplate"],
    args=sys.argv[2:]
  )

t.run()
