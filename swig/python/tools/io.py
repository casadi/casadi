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
      
      
@contextlib.contextmanager
def capture_stdout():
    import sys
    from cStringIO import StringIO
    oldout,olderr = sys.stdout, sys.stderr
    try:
        out=[StringIO(), StringIO()]
        sys.stdout.flush()
        sys.stderr.flush()
        sys.stdout,sys.stderr = out
        yield out
    finally:
        sys.stdout,sys.stderr = oldout, olderr
        out[0] = out[0].getvalue()
        out[1] = out[1].getvalue()
