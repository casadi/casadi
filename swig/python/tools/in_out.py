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
import contextlib

@contextlib.contextmanager
def nice_stdout():
  """
  
  There a various tricks in Python to capture redirect / capture sys.stdout.
  However, the C++ code is not aware of the Python sys.stdout
  
  This contextmanager captures the C++ output in memory and dumps it on the Python sys.stdout
  
  Caution:
    All C++ output is dumped to a pipe and only passed to stdout at the end of the call.
    This means that e.g. NLP iterates will not show up interactively.
    
    This could in theory be overcome by spawning a sister thread that periodically reads from the buffer and dumps to the Python stdout
    
    
  Usage:
  
  from casadi.tools import *

  x = SX.sym("x")

  with capture_stdout() as out:
    with nice_stdout():
      print "foo"
      x.sparsity().spy()
      
  """
  import os, sys
  (r,w) = os.pipe()
  try:
    fcntl.fcntl(r, fcntl.F_SETFL, os.O_NONBLOCK) # Make nonblocking, only on Linux
  except:
    pass
  sys.stdout.flush()
  backup = os.dup(1)
  os.dup2(w, 1)
  try:
      yield
  finally:
      os.dup2(backup, 1)
      os.write(w,"x")
      sys.stdout.write(os.read(r,2**20)[:-1])
      os.close(r)
      os.close(w)
      os.close(backup)
      
      
@contextlib.contextmanager
def capture_stdout():
    import sys
    try:
      from cStringIO import StringIO
    except:
      from io import StringIO
    oldout,olderr = sys.stdout, sys.stderr
    try:
        out=[StringIO(), StringIO()]
        sys.stdout,sys.stderr = out
        yield out
    finally:
        sys.stdout,sys.stderr = oldout, olderr
        out[0] = out[0].getvalue()
        out[1] = out[1].getvalue()
