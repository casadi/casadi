import sys
src = sys.argv[1]
import os
import shutil

from testsuite import TestSuite

def setdummybackend(dir):
  """
  Set the matplotlib backend to a dummy value, so that no graphs are generated.
  """
  shutil.copy('internal/pylab.py',dir)
  f=file(dir+'/matplotlibrc','w')
  f.write('backend : Template')
  f.close()

def removedummybackend(dir):
  if os.path.exists(dir+'/matplotlibrc'):
    os.remove(dir+'/matplotlibrc')
  if os.path.exists(dir+'/pylab.pyc'):
    os.remove(dir+'/pylab.pyc')
  if os.path.exists(dir+'/pylab.py'):
    os.remove(dir+'/pylab.py')
    
t = TestSuite(dirname=src,
  suffix="py",
  preRun=setdummybackend,
  postRun=removedummybackend,
  command = lambda dir,fn:  ["python", fn],
  skipdirs=[".svn","ctemplate"],
    args=sys.argv[2:]
  )
  
t.run()
